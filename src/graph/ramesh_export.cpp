#include "ramesh.h"
#include <fstream>
#include <iomanip>

namespace RaMesh {
    /* —— 供两个 export… 调用的内部工具 —— */
    bool emitMafBlock(std::ostream& os,
        const BlockPtr& blk,
        const std::map<SpeciesName,
        SeqPro::SharedManagerVariant>& seq_mgrs,
        bool  only_primary,
        bool  pairwise_mode,
        bool  allow_reverse /* = true */)
    {
        if (!blk) return false;

        /* ---------- 1. 收集 segment 列表 ---------- */
        struct Rec { SpeciesName sp; ChrName chr; SegPtr seg; };
        std::vector<Rec> recs;             recs.reserve(blk->anchors.size());
        bool have_reverse = false;

        for (auto& [sp_chr, seg] : blk->anchors) {
            if (only_primary && !seg->isPrimary()) continue;
            if (!allow_reverse && seg->strand == REVERSE) have_reverse = true;
            recs.push_back({ sp_chr.first, sp_chr.second, seg });
        }
        if (recs.size() < 2 || (!allow_reverse && have_reverse)) return false;

        /* ---------- 2. 参考行排第一 ---------- */
        auto it_ref = std::find_if(recs.begin(), recs.end(),
            [&](const Rec& r) { return r.chr == blk->ref_chr; });
        if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);

        /* ---------- 3. 准备 lambda：fetchSeq / fetchLen ---------- */
        auto fetchSeq = [](const SeqPro::ManagerVariant& mv,
            const ChrName& chr, Coord_t b, Coord_t l) {
                return std::visit([&](auto& p) { return p->getSubSequence(chr, b, l); }, mv);
            };
        auto fetchLen = [](const SeqPro::ManagerVariant& mv,
            const ChrName& chr) -> uint64_t {
                return std::visit([&](auto& p) { return p->getSequenceLength(chr); }, mv);
            };

        /* ---------- 4. 提取原始子串 ---------- */
        std::unordered_map<ChrName, std::string> seqs;
        std::unordered_map<ChrName, Cigar_t>     cigars;

        for (const Rec& r : recs) {
            auto mit = seq_mgrs.find(r.sp);
            if (mit == seq_mgrs.end()) return false;           // manager 缺失

            std::string raw = fetchSeq(*mit->second, r.chr,
                r.seg->start, r.seg->length);
            if (r.seg->strand == Strand::REVERSE) reverseComplement(raw);

            ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
            seqs.emplace(key, std::move(raw));
            cigars.emplace(key, r.seg->cigar);
        }

        /* ---------- 5. 归并对齐 ---------- */
        const ChrName ref_key = pairwise_mode ?
            recs.front().chr :
            recs.front().sp + "." + recs.front().chr;
        try {
            mergeAlignmentByRef(ref_key, seqs, cigars);
        }
        catch (const std::exception& e) {
            spdlog::warn("mergeAlignmentByRef failed: {}", e.what());
            return false;
        }

        /* ---------- 6. 写块 ---------- */
        os << "a score=0\n";

        auto write_row = [&](const Rec& r, const std::string& aln) {
            const auto& mgr = *seq_mgrs.at(r.sp);
            uint64_t chr_len = fetchLen(mgr, r.chr);

            os << "s "
                << std::left << std::setw(20)
                << (pairwise_mode ? r.chr : r.sp + "." + r.chr)
                << std::right << std::setw(12) << r.seg->start
                << std::setw(12) << r.seg->length
                << ' ' << (r.seg->strand == Strand::FORWARD ? '+' : '-')
                << std::setw(12) << chr_len
                << ' ' << aln << '\n';
            };

        write_row(recs.front(), seqs.at(ref_key));
        for (std::size_t i = 1; i < recs.size(); ++i) {
            const Rec& r = recs[i];
            ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
            write_row(r, seqs.at(key));
        }
        os << '\n';
    }

    /* 允许正反链 */
    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());

        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";

        for (const auto& wblk : blocks)
            emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode,
                /*allow_reverse=*/true);
    }

    /* 过滤掉含反向片段的块 */
    void RaMeshMultiGenomeGraph::exportToMafWithoutReverse(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool only_primary,
        bool pairwise_mode) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty()) fs::create_directories(maf_path.parent_path());

        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";

        for (const auto& wblk : blocks)
            emitMafBlock(ofs, wblk.lock(), seq_mgrs, only_primary, pairwise_mode,
                /*allow_reverse=*/false);
    }

    

} // namespace RaMesh

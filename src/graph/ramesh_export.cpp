#include "ramesh.h"
#include <fstream>
#include <iomanip>

namespace RaMesh {

    // TODO 目前还不支持cigar比对结果
    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::ManagerVariant>& seqpro_managers,
        bool                                                         only_primary,
        bool                                                         is_pairwise) const
    {
        namespace fs = std::filesystem;
        if (!maf_path.parent_path().empty())
            fs::create_directories(maf_path.parent_path());

        std::ofstream ofs(maf_path.string(),
            std::ios::out | std::ios::trunc | std::ios::binary);
        if (!ofs.is_open())
            throw std::runtime_error("Cannot open MAF file: " + maf_path.string());

        ofs << "##maf version=1 scoring=none\n";

        const auto revComp = [](std::string s) {
            auto comp = [](char c) -> char {
                switch (std::toupper(c)) {
                case 'A': return 'T'; case 'T': return 'A';
                case 'C': return 'G'; case 'G': return 'C';
                default: return 'N';
                }
                };
            for (char& c : s) c = comp(c);
            std::reverse(s.begin(), s.end());
            return s;
            };

        std::shared_lock gLock(rw);
        for (const auto& weak_blk : blocks)
        {
            BlockPtr blk = weak_blk.lock();
            if (!blk) continue;

            std::shared_lock bLock(blk->rw);

            /* ---------- 收集全部（物种, 染色体, Segment） ---------- */
            struct Rec { SpeciesName sp; ChrName chr; SegPtr seg; };
            std::vector<Rec> recs;
            recs.reserve(blk->anchors.size());

            for (const auto& [sp_chr, seg] : blk->anchors) {
                if (only_primary && !seg->isPrimary()) continue;
                recs.push_back({ sp_chr.first, sp_chr.second, seg });
            }
            if (recs.size() < 2) continue;     // 至少要有两条才能成块

            /* ---------- 让参考染色体行位于 recs[0] ---------- */
            auto it_ref = std::find_if(recs.begin(), recs.end(),
                [&](const Rec& r) { return r.chr == blk->ref_chr; });
            if (it_ref != recs.end() && it_ref != recs.begin())
                std::swap(*recs.begin(), *it_ref);           // 参考行放最前

            ofs << "a score=0\n";

            /* ---------- 写出所有 s-record ---------- */
            for (const auto& r : recs)
            {
                const SegPtr seg = r.seg;

                /* 1) 拿到对应序列管理器 */
                const auto mit = seqpro_managers.find(r.sp);
                if (mit == seqpro_managers.end())
                    continue;                                // 缺 manager：跳过当前记录

                uint64_t chr_len = 0;
                std::string subseq;

                std::visit([&](auto const& up) {
                    if (!up) return;
                    using PtrT = std::decay_t<decltype(up)>;
                    if constexpr (std::is_same_v<typename PtrT::element_type,
                        SeqPro::SequenceManager>)
                    {
                        chr_len = up->getSequenceLength(r.chr);
                        subseq = up->getSubSequence(r.chr, seg->start, seg->length);
                    }
                    else /* MaskedSequenceManager */
                    {
                        const auto& ori = up->getOriginalManager();
                        chr_len = ori.getSequenceLength(r.chr);
                        subseq = ori.getSubSequence(r.chr, seg->start, seg->length);
                    }
                    }, mit->second);

                if (subseq.empty()) continue;                // 无序列：跳过

                uint64_t maf_start = seg->start;
                if (seg->strand == Strand::REVERSE) {
                    maf_start = chr_len - (seg->start + seg->length);
                    reverseComplement(subseq);
                }

                ofs << "s "
                    << std::left << std::setw(20) << (is_pairwise ? r.chr : r.sp + "." + r.chr)
                    << std::right << std::setw(12) << maf_start
                    << std::setw(12) << seg->length
                    << ' ' << (seg->strand == Strand::FORWARD ? '+' : '-')
                    << std::setw(12) << chr_len
                    << ' ' << subseq
                    << '\n';
            }
            ofs << '\n';
        }
    }


    

} // namespace RaMesh

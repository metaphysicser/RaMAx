#include "ramesh.h"
#include <fstream>
#include <iomanip>

namespace RaMesh {

    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::ManagerVariant>& seq_mgrs,
        bool  only_primary,
        bool  pairwise_mode) const
    {
        namespace fs = std::filesystem;

        /* -------- 0. 打开文件、写头 -------- */
        if (!maf_path.parent_path().empty())
            fs::create_directories(maf_path.parent_path());

        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";

        /* -------- 小工具 -------- */
        auto fetchSeq = [&](const SeqPro::ManagerVariant& mv,
            const ChrName& chr, uint64_t b, uint64_t l)->std::string {
                return std::visit([&](auto& p) { return p->getSubSequence(chr, b, l);}, mv);
            };
        auto fetchLen = [&](const SeqPro::ManagerVariant& mv,
            const ChrName& chr)->uint64_t {
                return std::visit([&](auto& p) { return p->getSequenceLength(chr);}, mv);
            };

        /* ====================================================== */
        for (auto const& wblk : blocks)
        {
            BlockPtr blk = wblk.lock();
            if (!blk) continue;

            /* 1. 收集 Rec 列表 */
            struct Rec { SpeciesName sp; ChrName chr; SegPtr seg; };
            std::vector<Rec> recs;
            recs.reserve(blk->anchors.size());
            for (auto const& [sp_chr, seg] : blk->anchors) {
                if (only_primary && !seg->isPrimary()) continue;
                recs.push_back({ sp_chr.first, sp_chr.second, seg });
            }
            if (recs.size() < 2) continue;

            /* 2. 把参考行放 recs[0] */
            auto it_ref = std::find_if(recs.begin(), recs.end(),
                [&](const Rec& r) { return r.chr == blk->ref_chr; });
            if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);

            const Rec& refRec = recs.front();
            const SegPtr refSeg = refRec.seg;

            /* 3. 拿参考序列 & 长度 */
            auto itMgr = seq_mgrs.find(refRec.sp);
            if (itMgr == seq_mgrs.end()) continue;
            std::string ref_raw = fetchSeq(itMgr->second, refRec.chr, refSeg->start, refSeg->length);
            uint64_t    ref_chr_len = fetchLen(itMgr->second, refRec.chr);
            if (refSeg->strand == Strand::REVERSE) reverseComplement(ref_raw);

            /* 4. 写块头 */
            ofs << "a score=0\n";

            /* ===== 遍历 recs，先写参考行，再写每条 query 行 ===== */
            bool ref_written = false;
            for (const Rec& r : recs)
            {
                /* 4-a. Manager */
                auto mgrIt = seq_mgrs.find(r.sp);
                if (mgrIt == seq_mgrs.end()) continue;
                const SegPtr seg = r.seg;

                /* 4-b. 拿 query 序列 & 长度 */
                std::string qry_raw = fetchSeq(mgrIt->second, r.chr, seg->start, seg->length);
                uint64_t    qry_chr_len = fetchLen(mgrIt->second, r.chr);
                if (seg->strand == Strand::REVERSE) reverseComplement(qry_raw);

                /* 4-c. 生成对齐串 */
                auto [ref_aln, qry_aln] = buildAlignment(ref_raw, qry_raw, seg->cigar);
                if (ref_aln.empty() || qry_aln.empty()) continue;        // 保护

                /* 4-d. 如未写过参考行 → 写一次带 gap 的参考行 */
                if (!ref_written)
                {
                    uint64_t ref_maf_start = (refSeg->strand == Strand::FORWARD)
                        ? refSeg->start
                        : ref_chr_len - (refSeg->start + refSeg->length);

                    ofs << "s "
                        << std::left << std::setw(20)
                        << (pairwise_mode ? refRec.chr : refRec.sp + "." + refRec.chr)
                        << std::right << std::setw(12) << ref_maf_start
                        << std::setw(12) << refSeg->length
                        << ' ' << (refSeg->strand == Strand::FORWARD ? '+' : '-')
                        << std::setw(12) << ref_chr_len
                        << ' ' << ref_aln << '\n';

                    ref_written = true;
                }

                /* 4-e. 若当前就是参考记录 → 跳过写 query 行 */
                if (&r == &refRec) continue;

                /* 4-f. 写 query 行（带 gap） */
                uint64_t qry_maf_start = (seg->strand == Strand::FORWARD)
                    ? seg->start
                    : qry_chr_len - (seg->start + seg->length);

                ofs << "s "
                    << std::left << std::setw(20)
                    << (pairwise_mode ? r.chr : r.sp + "." + r.chr)
                    << std::right << std::setw(12) << qry_maf_start
                    << std::setw(12) << seg->length
                    << ' ' << (seg->strand == Strand::FORWARD ? '+' : '-')
                    << std::setw(12) << qry_chr_len
                    << ' ' << qry_aln << '\n';
            }
            ofs << '\n';     // block 分隔
        }
    }

    

} // namespace RaMesh

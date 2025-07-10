#include "ramesh.h"
#include <fstream>
#include <iomanip>

namespace RaMesh {


    void RaMeshMultiGenomeGraph::exportToMaf(
        const FilePath& maf_path,
        const std::map<SpeciesName, SeqPro::SharedManagerVariant>& seq_mgrs,
        bool  only_primary,
        bool  pairwise_mode) const
    {
        namespace fs = std::filesystem;

        /* ---- 0. 打开文件、写头 ---- */
        if (!maf_path.parent_path().empty())
            fs::create_directories(maf_path.parent_path());

        std::ofstream ofs(maf_path, std::ios::binary | std::ios::trunc);
        if (!ofs) throw std::runtime_error("Cannot open: " + maf_path.string());
        ofs << "##maf version=1 scoring=none\n";

        /* ---- 序列提取工具 ---- */
/* ---- 统一使用“原始管理器”取子串 ---- */
        auto fetchSeq = [&](const SeqPro::ManagerVariant& mv,
            const ChrName& chr, uint64_t b, uint64_t l) -> std::string
            {
                return std::visit([&](auto& up) -> std::string {
                    using PtrT = std::decay_t<decltype(up)>;

                    if constexpr (std::is_same_v<PtrT,
                        std::unique_ptr<SeqPro::SequenceManager>>) {
                        /* 直接 SequenceManager */
                        return up->getSubSequence(chr, b, l);

                    }
                    else { // → MaskedSequenceManager
                    //    /* 回退到底层原始管理器 */
                    return up->getOriginalManager()
                            .getSubSequence(chr, b, l);
                    }
                    }, mv);
            };

        /* ---- 同理：取染色体长度 ---- */
        auto fetchLen = [&](const SeqPro::ManagerVariant& mv,
            const ChrName& chr) -> uint64_t
            {
                return std::visit([&](auto& up) -> uint64_t {
                    using PtrT = std::decay_t<decltype(up)>;

                    if constexpr (std::is_same_v<PtrT,
                        std::unique_ptr<SeqPro::SequenceManager>>) {
                        return up->getSequenceLength(chr);
                    }
                    else {
                    return up->getOriginalManager()
                            .getSequenceLength(chr);
                    }
                    return up->getSequenceLength(chr);
                    }, mv);

            };


        /* ========================================================= */
        for (const auto& wblk : blocks)
        {
            BlockPtr blk = wblk.lock();
            if (!blk) continue;

            //---------------- 1. 收集 segment 列表 -----------------
            struct Rec { SpeciesName sp; ChrName chr; SegPtr seg; };
            std::vector<Rec> recs;
            recs.reserve(blk->anchors.size());

            for (auto& [sp_chr, seg] : blk->anchors) {
                if (only_primary && !seg->isPrimary()) continue;
                recs.push_back({ sp_chr.first, sp_chr.second, seg });
            }
            if (recs.size() < 2) continue;

            //---------------- 2. 让参考排第一个 --------------------
            auto it_ref = std::find_if(recs.begin(), recs.end(),
                [&](const Rec& r) { return r.chr == blk->ref_chr; });
            if (it_ref != recs.end()) std::swap(*recs.begin(), *it_ref);

            const Rec& refRec = recs.front();
            const SegPtr refSeg = refRec.seg;

            //---------------- 3. 提取所有原始子串 -------------------
            std::unordered_map<ChrName, std::string> seqs;   // key -> rawSeq
            std::unordered_map<ChrName, Cigar_t>     cigars; // key -> CIGAR (query 侧)

            for (const Rec& r : recs)
            {
                /* 3-a 找到序列管理器 */
                auto mgrIt = seq_mgrs.find(r.sp);
                if (mgrIt == seq_mgrs.end()) { seqs.clear(); break; }

                /* 3-b 获取并按 strand 处理序列 */
                const SegPtr seg = r.seg;

                std::string raw = fetchSeq(*mgrIt->second, r.chr, seg->start, seg->length);
                if (seg->strand == Strand::REVERSE) reverseComplement(raw);

                /* 3-c key 名称（也是 MAF 行名称） */
                ChrName key_name = pairwise_mode ? r.chr : r.sp + "." + r.chr;
                seqs.emplace(key_name, std::move(raw));

                
                cigars.emplace(key_name, seg->cigar);
            }
            if (seqs.empty()) continue;           // 缺管理器则跳过

            //---------------- 4. 归并成多序列对齐 ------------------
            const ChrName ref_key = pairwise_mode ? refRec.chr : refRec.sp + "." + refRec.chr;
            if (refRec.sp != "simChimp") {
                std::cout << "";
            }
            // try {
            //     mergeAlignmentByRef(ref_key, seqs, cigars);     // 就地修改 seqs
            // }
            // catch (const std::exception& e) {
            //     spdlog::warn("mergeAlignmentByRef failed: {}", e.what());
            //     continue;
            // }

            //---------------- 5. 写 MAF 块头 ----------------------
            ofs << "a score=0\n";

            /* 方便后面统一写行的 lambda */
            auto write_row = [&](const Rec& r, const std::string& aln) {
                auto mgrIt = seq_mgrs.find(r.sp);
                uint64_t chr_len = fetchLen(*mgrIt->second, r.chr);


                ofs << "s "
                    << std::left << std::setw(20)
                    << (pairwise_mode ? r.chr : r.sp + "." + r.chr)
                    << std::right << std::setw(12) << r.seg->start
                    << std::setw(12) << r.seg->length            // 非 gap 长度
                    << ' ' << (r.seg->strand == Strand::FORWARD ? '+' : '-')
                    << std::setw(12) << chr_len
                    << ' ' << aln << '\n';
                };

            //---------------- 6. 写参考行 -------------------------
            write_row(refRec, seqs[ref_key]);

            //---------------- 7. 写 query 行 ----------------------
            for (std::size_t i = 1; i < recs.size(); ++i) {
                const Rec& r = recs[i];
                ChrName key = pairwise_mode ? r.chr : r.sp + "." + r.chr;
                write_row(r, seqs[key]);
            }

            ofs << '\n';   // 块间空行
        }
    }


    

} // namespace RaMesh

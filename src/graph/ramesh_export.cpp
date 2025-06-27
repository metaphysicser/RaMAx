#include "ramesh.h"
#include <fstream>
#include <iomanip>

namespace RaMesh {

    void RaMeshMultiGenomeGraph::exportToMaf(const FilePath& maf_path,
        bool only_primary) const
    {
		spdlog::info("Exporting MAF to: {}", maf_path.string());
        auto dir = maf_path.parent_path();

        if (!dir.empty() && !std::filesystem::exists(dir)) {
            if (!std::filesystem::create_directories(dir)) {
                throw std::runtime_error("Failed to create dir: " + dir.string());
            }
        }

        // —— 用 .string() 显式打开 —— 
        std::ofstream ofs(maf_path.string(),
            std::ios::out | std::ios::trunc | std::ios::binary);

        if (!ofs) {
            throw std::runtime_error("Cannot open MAF file: " + maf_path.string());
        }

        // —— header —— 
        ofs << "##maf version=1 scoring=none\n";
        ofs.flush();  // 确保 header 刷出到磁盘

        std::shared_lock gLock(rw);          // 多线程只读即可
        for (const auto& weak_blk : blocks)
        {
            BlockPtr blk = weak_blk.lock();
            if (!blk) continue;

            std::shared_lock bLock(blk->rw);

            // —— 只导出二元块（简化示例）——
            if (blk->anchors.size() != 2) continue;

            struct Rec { SpeciesName sp; ChrName chr; SegPtr seg; };
            std::vector<Rec> recs;
            for (const auto& [sp_chr, seg] : blk->anchors) {
                if (only_primary && !seg->isPrimary()) continue;
                recs.push_back({ sp_chr.first, sp_chr.second, seg });
            }
            if (recs.size() != 2) continue;

            // —— 输出一个 MAF block —— 
            ofs << "a score=0\n";
            if (!std::filesystem::exists(maf_path.string()) || std::filesystem::file_size(maf_path.string()) == 0) {
                throw std::runtime_error("Header not written – path="
                    + maf_path.string());
            }
            for (const auto& r : recs)
            {
                const SegPtr seg = r.seg;
                const uint64_t chr_len = 0;              // TODO: 若能获取真实长度请替换
                uint64_t maf_start = seg->start;         // forward strand坐标
                if (seg->strand == Strand::REVERSE && chr_len)
                    maf_start = chr_len - (seg->start + seg->length);

                ofs << "s "
                    << std::left << std::setw(20) << (r.sp + "." + r.chr)  // src
                    << std::right << std::setw(12) << maf_start             // start
                    << std::setw(12) << seg->length                         // size
                    << ' ' << (seg->strand == Strand::FORWARD ? '+' : '-')  // strand
                    << std::setw(12) << chr_len                             // srcSize
                    << ' ' << std::string(seg->length, 'N')                 // seq (占位)
                    << '\n';
            }
            ofs << '\n';
        }
        ofs.close();
    }
    

} // namespace RaMesh

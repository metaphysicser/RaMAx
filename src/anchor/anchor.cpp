#include "anchor.h"
#include "data_process.h" 

// 构造函数（注释掉了）：初始化参考物种、查询物种及其锚点集合，并重建 R 树。
//PairGenomeAnchor::PairGenomeAnchor(SpeciesName ref,
//    SpeciesName query,
//    AnchorPtrListVec init_anchors)
//    : ref_species(ref)
//    , query_species(query)
//    , anchors(std::move(init_anchors))
//{
//    rebuild(); // 构建 R 树索引
//}

// 重建 R 树索引：将 anchors 中所有 Anchor 加入 R 树中。
//
//void PairGenomeAnchor::rebuild() {
//    // 1) 清空已有的 R 树索引
//    anchor_rtree.clear();
//
//    // 2) 遍历 anchors 的每个 list，插入 AnchorPtr 到 R-tree 中
//    for (auto const& list_group : anchors) {
//        for (auto const& ap : list_group) {
//            anchor_rtree.insert({ make_box(*ap), ap });
//        }
//    }
//}

// 将一个 Anchor 插入 R 树中
//void PairGenomeAnchor::insertAnchorRtree(AnchorPtr ap) {
//    anchor_rtree.insert({ make_box(*ap), ap });
//}

// 从 R 树中移除一个 Anchor
//void PairGenomeAnchor::removeAnchorRtree(AnchorPtr ap) {
//    // 从 R 树中删除 Anchor
//    anchor_rtree.remove({ make_box(*ap), ap });
//}

// 查询与给定区域相交的所有锚点
//std::vector<AnchorPtr> PairGenomeAnchor::query(const AnchorBox& region) const {
//    // 在 R 树中查找与 region 相交的 AnchorEntry
//    std::vector<AnchorEntry> hits;
//    anchor_rtree.query(bgi::intersects(region),
//        std::back_inserter(hits));
//
//    // 提取命中的 AnchorPtr 到输出向量中
//    std::vector<AnchorPtr> out;
//    out.reserve(hits.size());
//    for (auto const& e : hits) {
//        out.push_back(e.second);
//    }
//    return out;
//}

// ------------------------------------------------------------------
// 保存锚点数据到文件：将二维列表结构 anchors 写入二进制文件
// ------------------------------------------------------------------
bool saveAnchors(const std::string& filename,
    const AnchorPtrListVec& anchors)
{
    std::ofstream os(filename, std::ios::binary); // 以二进制方式打开文件
    if (!os) return false; // 打开失败

    cereal::BinaryOutputArchive oar(os); // 创建 cereal 的输出归档对象

    // 1) 写入外层 vector 的大小
    uint64_t outer = anchors.size();
    oar(outer);

    // 2) 对每个内部 list：
    for (const auto& lst : anchors) {
        uint64_t inner = lst.size(); // 写入 list 大小
        oar(inner);

        // 写入每个 Anchor 对象的内容
        for (const auto& ap : lst) {
            oar(*ap); // 保存对象内容，而不是指针本身
        }
    }
    return static_cast<bool>(os); // 返回写入是否成功
}

// ------------------------------------------------------------------
// 从文件读取锚点数据：与 saveAnchors 格式对应
// ------------------------------------------------------------------
bool loadAnchors(const std::string& filename,
    AnchorPtrListVec& anchors)
{
    std::ifstream is(filename, std::ios::binary); // 打开二进制文件
    if (!is) return false;

    cereal::BinaryInputArchive iar(is); // 创建 cereal 的输入归档对象

    uint64_t outer;
    iar(outer); // 读取外层 vector 大小

    anchors.clear();
    anchors.reserve(static_cast<size_t>(outer)); // 预分配空间

    // 读取每个内部 list
    for (uint64_t i = 0; i < outer; ++i) {
        uint64_t inner;
        iar(inner); // 读取内部 list 大小

        AnchorPtrList lst;
        for (uint64_t j = 0; j < inner; ++j) {
            auto ap = std::make_shared<Anchor>(); // 创建新的 Anchor 智能指针
            iar(*ap); // 从文件读取 Anchor 内容到 ap 中
            lst.push_back(std::move(ap));
        }
        anchors.emplace_back(std::move(lst));
    }
    return static_cast<bool>(is); // 返回读取是否成功
}

// ------------------------------------------------------------------
// (1) 保存 std::vector<AnchorPtrListVec> 到文件
// ------------------------------------------------------------------
bool saveAnchorsSets(const std::string& filename,
    const std::vector<AnchorPtrListVec>& all_sets)
{
    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    /* ---------- 写入最外层大小 ---------- */
    uint64_t lvl0 = all_sets.size();
    oar(lvl0);

    /* ---------- 逐层写入 ---------- */
    for (const auto& anchors : all_sets)
    {
        /* ---- 第二层：AnchorPtrListVec 大小 ---- */
        uint64_t lvl1 = anchors.size();
        oar(lvl1);

        for (const auto& lst : anchors)
        {
            /* -- 第三层：AnchorPtrList 大小 -- */
            uint64_t lvl2 = lst.size();
            oar(lvl2);

            for (const auto& ap : lst)
                oar(*ap);          // 仅写 Anchor 对象内容
        }
    }
    return static_cast<bool>(os);
}

// ------------------------------------------------------------------
// (2) 从文件读取 std::vector<AnchorPtrListVec>
// ------------------------------------------------------------------
bool loadAnchorsSets(const std::string& filename,
    std::vector<AnchorPtrListVec>& all_sets)
{
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t lvl0;
    iar(lvl0);                         // 读取最外层大小

    all_sets.clear();
    all_sets.reserve(static_cast<size_t>(lvl0));

    for (uint64_t i = 0; i < lvl0; ++i)
    {
        uint64_t lvl1;
        iar(lvl1);                     // 读取第二层大小

        AnchorPtrListVec anchors;
        anchors.reserve(static_cast<size_t>(lvl1));

        for (uint64_t j = 0; j < lvl1; ++j)
        {
            uint64_t lvl2;
            iar(lvl2);                 // 读取第三层大小

            AnchorPtrList lst;
            for (uint64_t k = 0; k < lvl2; ++k)
            {
                auto ap = std::make_shared<Anchor>();
                iar(*ap);              // 反序列化到对象
                lst.push_back(std::move(ap));
            }
            anchors.emplace_back(std::move(lst));
        }
        all_sets.emplace_back(std::move(anchors));
    }
    return static_cast<bool>(is);
}

bool saveAnchorVec3D(const std::string& filename, const AnchorVec3DPtr& data) {
    if (!data) return false;

    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    uint64_t dim0 = data->size();  // 最外层维度
    oar(dim0);

    for (const auto& layer : *data) {
        uint64_t dim1 = layer.size();  // 中间层维度
        oar(dim1);

        for (const auto& vec : layer) {
            uint64_t dim2 = vec.size();  // 最内层 Anchor 数
            oar(dim2);

            for (const auto& anchor : vec) {
                oar(anchor);  // 注意不是指针，这里是值类型 Anchor
            }
        }
    }

    return static_cast<bool>(os);
}

bool loadAnchorVec3D(const std::string& filename, AnchorVec3DPtr& data) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t dim0;
    iar(dim0);

    data = std::make_shared<AnchorVec3D>();
    data->resize(dim0);

    for (uint64_t i = 0; i < dim0; ++i) {
        uint64_t dim1;
        iar(dim1);

        (*data)[i].resize(dim1);

        for (uint64_t j = 0; j < dim1; ++j) {
            uint64_t dim2;
            iar(dim2);

            AnchorVec vec;
            vec.reserve(dim2);

            for (uint64_t k = 0; k < dim2; ++k) {
                Anchor a;
                iar(a);
                vec.push_back(std::move(a));
            }

            (*data)[i][j] = std::move(vec);
        }
    }

    return static_cast<bool>(is);
}


//void groupMatchByQueryRef(MatchVec3DPtr& anchors,
//    MatchByQueryRef& unique_anchors,
//    MatchByQueryRef& repeat_anchors,
//    FastaManager& ref_fasta_manager,
//    FastaManager& query_fasta_manager,
//    uint_t thread_num) {
//    //----------------------------------------------------------------
//    // 0) 初始化二维目标矩阵
//    //----------------------------------------------------------------
//    const uint_t ref_chr_cnt = ref_fasta_manager.idx_map.size();
//    const uint_t query_chr_cnt = query_fasta_manager.idx_map.size();
//
//    // resize 外层：query 维度
//    unique_anchors.resize(query_chr_cnt);
//    repeat_anchors.resize(query_chr_cnt);
//
//    // resize 内层：ref 维度
//    for (uint_t q = 0; q < query_chr_cnt; ++q) {
//        unique_anchors[q].resize(ref_chr_cnt);   // AnchorsByRef
//        repeat_anchors[q].resize(ref_chr_cnt);   // AnchorsByRef
//    }
//
//    //----------------------------------------------------------------
//    // 1) 行级互斥锁：每个 query-chr 一把
//    //----------------------------------------------------------------
//    std::vector<std::mutex> rowLocks(query_chr_cnt);
//
//    //----------------------------------------------------------------
//    // 2) 并行遍历 3D 数据，每个 slice 启动一个任务
//    //----------------------------------------------------------------
//    ThreadPool pool(thread_num);
//
//    for (auto& slice : *anchors)              // slice: std::vector<AnchorVec>
//    {
//        pool.enqueue([&, &slice]() {
//
//            if (slice.empty()) return;
//
//            // 2a) 获取 slice 的 query-chr 索引（所有 vec 同一 chr）
//            const Match& first = slice.front().front();
//            auto itQ = query_fasta_manager.idx_map.find(
//                first.query_region.chr_name);
//            if (itQ == query_fasta_manager.idx_map.end()) return;
//            const uint_t qIdx = itQ->second;
//
//            // 2b) 遍历 slice 中每个 AnchorVec
//            for (auto& vec : slice) {
//                if (vec.empty()) continue;
//
//                const ChrName& rName = vec.front().ref_region.chr_name;
//                auto itR = ref_fasta_manager.idx_map.find(rName);
//                if (itR == ref_fasta_manager.idx_map.end()) continue;
//                const uint_t rIdx = itR->second;
//
//                /* ---- 临界区：只锁当前 query 行 ---- */
//                {
//                    std::lock_guard<std::mutex> lk(rowLocks[qIdx]);
//                    MatchVec& target =
//                        (vec.size() == 1)
//                        ? unique_anchors[qIdx][rIdx]
//                        : repeat_anchors[qIdx][rIdx];
//
//                        target.insert(target.end(),
//                            std::make_move_iterator(vec.begin()),
//                            std::make_move_iterator(vec.end()));
//                }
//            }
//            });
//    }
//
//    pool.waitAllTasksDone();
//
//    //----------------------------------------------------------------
//    // 3) 释放原始 3D 数据以节省内存
//    //----------------------------------------------------------------
//    anchors->clear();
//    anchors->shrink_to_fit();
//}

//--------------------------------------------------------------------
// 适配 4‑D 结构的聚类分发
//--------------------------------------------------------------------
void groupMatchByQueryRef(MatchVec3DPtr& anchors,
    MatchByStrandByQueryRef& unique_anchors,
    MatchByStrandByQueryRef& repeat_anchors,
    FastaManager& ref_fasta_manager,
    FastaManager& query_fasta_manager,
    uint_t thread_num)
{
    //------------------------------------------------------------
    // 0) 初始化 4‑D 输出矩阵  [strand][query][ref]
    //------------------------------------------------------------
    constexpr uint_t STRAND_CNT = 2;            // FORWARD=0, REVERSE=1

    const uint_t ref_chr_cnt = static_cast<uint_t>(ref_fasta_manager.idx_map.size());
    const uint_t query_chr_cnt = static_cast<uint_t>(query_fasta_manager.idx_map.size());

    auto initTarget = [&](MatchByStrandByQueryRef& tgt) {
        tgt.resize(STRAND_CNT);
        for (uint_t s = 0; s < STRAND_CNT; ++s) {
            tgt[s].resize(query_chr_cnt);              // [query]
            for (uint_t q = 0; q < query_chr_cnt; ++q)
                tgt[s][q].resize(ref_chr_cnt);         // [ref]
        }
        };

    initTarget(unique_anchors);
    initTarget(repeat_anchors);

    //------------------------------------------------------------
    // 1) 行级互斥锁: 针对 (strand, query) 组合各一把锁
    //------------------------------------------------------------
    std::vector<std::mutex> rowLocks(STRAND_CNT * query_chr_cnt);
    auto rowLockIdx = [&](uint_t strand, uint_t qIdx) {
        return strand * query_chr_cnt + qIdx;
        };

    //------------------------------------------------------------
    // 2) 并行遍历 3‑D 数据，slice 级别开任务
    //------------------------------------------------------------
    ThreadPool pool(thread_num);

    for (auto& slice : *anchors)
    {
        pool.enqueue([&, slice]() {              // ← 值捕获
            if (slice.empty()) return;

            const Match& first = slice.front().front();

            // 1. strand 索引归一化
            uint_t sIdx = (first.strand == REVERSE ? 1u : 0u);

            auto itQ = query_fasta_manager.idx_map.find(first.query_region.chr_name);
            if (itQ == query_fasta_manager.idx_map.end()) return;
            uint_t qIdx = itQ->second;

            for (auto& vec : slice)
            {
                if (vec.empty()) continue;
                auto itR = ref_fasta_manager.idx_map.find(vec.front().ref_region.chr_name);
                if (itR == ref_fasta_manager.idx_map.end()) continue;
                uint_t rIdx = itR->second;

                std::lock_guard<std::mutex> lk(rowLocks[rowLockIdx(sIdx, qIdx)]);
                MatchVec& tgt = (vec.size() == 1)
                    ? unique_anchors[sIdx][qIdx][rIdx]
                    : repeat_anchors[sIdx][qIdx][rIdx];
                    tgt.insert(tgt.end(),
                        std::make_move_iterator(vec.begin()),
                        std::make_move_iterator(vec.end()));
            }
            });
    }


    pool.waitAllTasksDone();

    //------------------------------------------------------------
    // 3) 释放原始 3‑D 数据以节省内存
    //------------------------------------------------------------
    anchors->clear();
    anchors->shrink_to_fit();
}

void sortMatchByRefStart(MatchByQueryRef& anchors, uint_t thread_num) {
    ThreadPool pool(thread_num);

    for (auto& row : anchors) {
        for (auto& vec : row) {
            pool.enqueue([&vec]() {
                if (vec.size() == 0) return;
                std::sort(vec.begin(), vec.end(), [](const Match& a, const Match& b) {
                    return (a.ref_region.start < b.ref_region.start) || (a.ref_region.start == b.ref_region.start && a.query_region.start < b.query_region.start);
                    });
                });
        }
    }

    pool.waitAllTasksDone();
}

void sortMatchByQueryStart(MatchByQueryRef& anchors, uint_t thread_num) {
    ThreadPool pool(thread_num);

    for (auto& row : anchors) {
        for (auto& vec : row) {
            pool.enqueue([&vec]() {
                if (vec.size() == 0) return;
                std::sort(vec.begin(), vec.end(), [](const Match& a, const Match& b) {
                    return a.query_region.start < b.query_region.start;
                    });
                });
        }
    }

    pool.waitAllTasksDone();
}

AnchorPtrVec findNonOverlapAnchors(const AnchorVec& anchors)
{
    AnchorPtrVec uniques;
    const size_t n = anchors.size();
    if (n == 0) return uniques;

    Coord_t max_end_so_far = 0;  // 记录到当前位置前的最大 end

    for (size_t i = 0; i < n; ++i) {
        const Anchor& curr = anchors[i];
        Coord_t curr_start = curr.match.ref_region.start;
        Coord_t curr_end = curr_start + curr.match.ref_region.length;

        bool overlap_left = (i > 0 && max_end_so_far > curr_start);
        bool overlap_right = (i + 1 < n &&
            anchors[i + 1].match.ref_region.start < curr_end);

        if (!overlap_left && !overlap_right)
            uniques.emplace_back(std::make_shared<Anchor>(curr));

        // 更新左侧最大 end
        if (curr_end > max_end_so_far) max_end_so_far = curr_end;
    }
    return std::move(uniques);
}

bool saveMatchVec3D(const std::string& filename, const MatchVec3DPtr& data) {
    if (!data) return false;

    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    uint64_t dim0 = data->size();  // 最外层
    oar(dim0);

    for (const auto& layer : *data) {
        uint64_t dim1 = layer.size();  // 第二层
        oar(dim1);

        for (const auto& vec : layer) {
            uint64_t dim2 = vec.size();  // 最内层
            oar(dim2);

            for (const auto& match : vec) {
                oar(match);  // 直接序列化 Match
            }
        }
    }

    return static_cast<bool>(os);
}

bool loadMatchVec3D(const std::string& filename, MatchVec3DPtr& data) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t dim0;
    iar(dim0);

    data = std::make_shared<MatchVec3D>();
    data->resize(dim0);

    for (uint64_t i = 0; i < dim0; ++i) {
        uint64_t dim1;
        iar(dim1);

        (*data)[i].resize(dim1);

        for (uint64_t j = 0; j < dim1; ++j) {
            uint64_t dim2;
            iar(dim2);

            MatchVec vec;
            vec.reserve(dim2);

            for (uint64_t k = 0; k < dim2; ++k) {
                Match m;
                iar(m);
                vec.push_back(std::move(m));
            }

            (*data)[i][j] = std::move(vec);
        }
    }

    return static_cast<bool>(is);
}





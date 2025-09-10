#include "align.h"

KSW2AlignConfig makeDefaultKSW2Config() {
    static int8_t simple_dna_mat[25];
    static bool initialized = false;
    if (!initialized) {
        // A C G T N -> 0 1 2 3 4
        const int match = 2;
        const int mismatch = -3;
        const int ambiguous = -1;  // 对N的惩罚较小

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                if (i == 4 || j == 4) {
                    simple_dna_mat[i * 5 + j] = ambiguous;  // N匹配
                }
                else if (i == j) {
                    simple_dna_mat[i * 5 + j] = match;
                }
                else {
                    simple_dna_mat[i * 5 + j] = mismatch;
                }
            }
        }
        initialized = true;
    }

    return {
        .mat = simple_dna_mat,
        .alphabet_size = 5,
        .gap_open = 8,
        .gap_extend = 1,
        .end_bonus = 0,
        .zdrop = 100,           // 较远匹配终止
        .band_width = -1,      // 合理的band可提高性能
        .flag = KSW_EZ_GENERIC_SC | KSW_EZ_RIGHT
    };
}

Cigar_t globalAlignKSW2(const std::string& ref,
    const std::string& query)
{  
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    // KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整
    KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};

    ksw_extz2_sse(0,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);


    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);           // KSW2 用 malloc()
    return cigar;
}

/**********************************************************************
*  extendAlignKSW2  ——  ends-free（seed-and-extend）比对
*    @param ref        参考片段（目标方向）
*    @param query      查询片段（同方向；若反链请先反向互补）
*    @param zdrop      Z-drop 剪枝阈值（默认 200）
*    @param band       带宽限制；<0 表示不限制
*    @return           Cigar_t（BAM 编码）
**********************************************************************/
Cigar_t extendAlignKSW2(const std::string& ref,
    const std::string& query,
    int zdrop)
{
    /* ---------- 1. 序列编码 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());
    for (size_t i = 0; i < ref.size(); ++i) ref_enc[i] = ScoreChar2Idx[(uint8_t)ref[i]];
    for (size_t i = 0; i < query.size(); ++i) qry_enc[i] = ScoreChar2Idx[(uint8_t)query[i]];

    /* ---------- 2. 配置 ---------- */
    KSW2AlignConfig cfg = makeTurboKSW2Config(query.size(), ref.size());
    cfg.zdrop = zdrop;       // 用于提前终止
    cfg.flag = KSW_EZ_EXTZ_ONLY     // ends-free extension
        | KSW_EZ_APPROX_MAX    // 跟踪 ez.max_q/max_t
        | KSW_EZ_APPROX_DROP   // 在 approximate 模式下触发 z-drop 就中断
        | KSW_EZ_RIGHT;        // （可选）gap 右对齐     // **关键**：启用 extension/ends-free
    // 若需要右对齐 gaps 建议保留 KSW_EZ_RIGHT
    cfg.end_bonus = 100;
    cfg.band_width = 64;

    /* ---------- 3. 调用 KSW2 ---------- */
    ksw_extz_t ez{};
    ksw_extz2_sse(nullptr,
        static_cast<int>(qry_enc.size()), qry_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        cfg.alphabet_size, cfg.mat,
        cfg.gap_open, cfg.gap_extend,
        cfg.band_width, cfg.zdrop, cfg.end_bonus,
        cfg.flag, &ez);

	// 赋值bool& if_zdrop,int& ref_end,int& qry_end
    /* ---------- 4. 拷贝 & 释放 ---------- */
    Cigar_t cigar;
    cigar.reserve(ez.n_cigar);
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    free(ez.cigar);                    // ksw2 使用 malloc
    return cigar;                      // 返回的 CIGAR 即延伸片段
}

// 依赖同一份 Cigar_t/type 定义、以及 WFA2 头文件
Cigar_t extendAlignWFA2(const std::string& ref,
    const std::string& query,
    int zdrop)        // Z-drop 阈值
{
    ///* ---------- 1. 配置 WFA 属性 ---------- */
    //wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;

    //// 打分模式：Affine gap
    //attr.distance_metric = gap_affine;
    //attr.affine_penalties.mismatch = 2;  // X
    //attr.affine_penalties.gap_opening = 3;  // O
    //attr.affine_penalties.gap_extension = 1;  // E

    //// 端自由：左端锚定，右端自由（右向延伸）
    //attr.alignment_form.span = alignment_endsfree;
    //attr.alignment_form.pattern_begin_free = 0;
    //attr.alignment_form.pattern_end_free = pattern_end_free;
    //attr.alignment_form.text_begin_free = 0;
    //attr.alignment_form.text_end_free = text_end_free;


    //// 内存/带宽启发式
    //attr.memory_mode = wavefront_memory_high;
    //attr.heuristic.strategy = wf_heuristic_zdrop;
    //attr.heuristic.zdrop = zdrop;   // 传入阈值
    //attr.heuristic.steps_between_cutoffs = 50;      // 每 50 波前检查一次
    //// 可选：自适应带宽剪枝（与全局版本一致）
    //attr.heuristic.strategy_secondary = wf_heuristic_banded_adaptive;
    //attr.heuristic.min_k = -10;
    //attr.heuristic.max_k = +10;

    ///* ---------- 2. 创建对齐器并执行 ---------- */
    //wavefront_aligner_t* wf_aligner = wavefront_aligner_new(&attr);
    //wavefront_align(wf_aligner,
    //    ref.c_str(), ref.length(),
    //    query.c_str(), query.length());

    ///* ---------- 3. 取出 CIGAR ---------- */
    //uint32_t* cigar_buf = nullptr;
    //int       cigar_len = 0;
    //cigar_get_CIGAR(wf_aligner->cigar,
    //    /*print=*/false,
    //    &cigar_buf,
    //    &cigar_len);

    ///* ---------- 4. 拷贝到用户侧结构 ---------- */
    //Cigar_t cigar;
    //cigar.reserve(cigar_len);
    //for (int i = 0; i < cigar_len; ++i)
    //    cigar.push_back(cigar_buf[i]);

    ///* ---------- 5. 清理 ---------- */
    //wavefront_aligner_delete(wf_aligner);
    //free(cigar_buf);   // `cigar_get_CIGAR` 内部 malloc

    //return cigar;      // 返回从序列左端到 “Z-drop / 序列右端” 的比对
}



Cigar_t globalAlignWFA2(const std::string& ref,
    const std::string& query)
{
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.mismatch = 2;      // X > 0
    attributes.affine_penalties.gap_opening = 3;   // O >= 0
    attributes.affine_penalties.gap_extension = 1; // E > 0
    attributes.memory_mode = wavefront_memory_high;
    attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
    attributes.heuristic.min_k = -10;
    attributes.heuristic.max_k = +10;
    attributes.heuristic.steps_between_cutoffs = 1;
    //// Create a WFAligner
    //
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);

    wavefront_align(wf_aligner, ref.c_str(), ref.length(), query.c_str(), query.length());
    /*wfa::WFAlignerGapAffine aligner(2, 3, 1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);

    aligner.alignEnd2End(ref, query);*/
    uint32_t* cigar_buffer; // Buffer to hold the resulting CIGAR operations.
    int cigar_length = 0; // Length of the CIGAR string.
    // Retrieve the CIGAR string from the wavefront aligner.
    // cigar_get_CIGAR(aligner->cigar, true, &cigar_buffer, &cigar_length);

    /* ---------- 4. 拷贝 / 释放 CIGAR ---------- */
    Cigar_t cigar;

    for (int i = 0; i < cigar_length; ++i)
        cigar.push_back(cigar_buffer[i]);

    wavefront_aligner_delete(wf_aligner);

    return cigar;
}

/* ──────────── 合并成 MSA (就地修改 seqs) ──────────── */
uint_t mergeAlignmentByRef(
    ChrName ref_name,
    std::unordered_map<ChrName, std::string>& seqs,
    const std::unordered_map<ChrName, Cigar_t>& cigars)
{
    auto ref_it = seqs.find(ref_name);
    if (ref_it == seqs.end())
        throw std::invalid_argument("mergeAlignmentByRef: ref not found");

    std::string& ref_raw = ref_it->second;
    uint_t total_aligned_length = ref_raw.size();
	RefAlignInfo insert_info;

    for (const auto& [key, cigar] : cigars) {
        if (key == ref_name) continue;
        auto q_it = seqs.find(key);
        if (q_it == seqs.end()) {
            throw std::invalid_argument("mergeAlignmentByRef: seq missing");
        }

		std::string& qry_raw = q_it->second;

		uint_t ref_pos = 0;
		uint_t qry_pos = 0;
        for (auto& unit : cigar) {
            uint32_t len;
            char op;
			intToCigar(unit, op, len);

            if (op == 'D') {
				qry_raw.insert(qry_pos, len, '-');
                ref_pos += len;
            }
            else if (op == 'I') {
                std::string ins = qry_raw.substr(qry_pos, len);

                // 2) 在 insert_info 里插入或更新
                auto it = insert_info.find(ref_pos);
                if (it != insert_info.end()) {
                    it->second.seqs[key] = ins;
                }
                else {
                    InsertInfo info;
                    info.seqs[key] = ins;
                    insert_info[ref_pos] = std::move(info);
                }

                // 3) 从原始 query 序列里移除这段已“消费”的子串
                qry_raw.erase(qry_pos, len);
            }
            else {
                ref_pos += len;
                qry_pos += len;
            }
        }
    }

    uint_t offset = 0;
	for (auto& [ref_pos, info] : insert_info) {
		info.alignSeqs(); // 对齐所有插入序列
		if (info.ref_name.empty()) continue; // 没有参考序列，跳过
		for (auto& [sp_name, seq] : seqs) {
			auto it = info.seqs.find(sp_name);
            if (it != info.seqs.end()) {
				seq.insert(ref_pos, it->second); // 在 ref_pos 位置插入
            }
            else {
                seq.insert(ref_pos, info.total_length, '-');   // 直接用 string::insert 重载
            }
		}
        total_aligned_length += info.total_length; // 更新总长度
	}

    for (auto& [chr, seq] : seqs) {
        if (seq.size() != total_aligned_length) {
            std::cout << "";
        }
    }

    return total_aligned_length;

}





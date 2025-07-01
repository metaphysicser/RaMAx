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





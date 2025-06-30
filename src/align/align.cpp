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
    const std::string& query,
    const KSW2AlignConfig& cfg_in)
{
    /* ---------- 1. 编码序列 ---------- */
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> qry_enc(query.size());

    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        qry_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    /* ---------- 2. 复制 cfg 并修正常见坑 ---------- */
    KSW2AlignConfig cfg = cfg_in;                   // 本地副本可调整

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



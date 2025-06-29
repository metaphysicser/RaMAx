#include "align.h"

KSW2AlignConfig makeDefaultKSW2Config() {
    static int8_t hoxd70_flat[25];
    static bool initialized = false;
    if (!initialized) {
        for (int i = 0; i < 5; ++i)
            for (int j = 0; j < 5; ++j)
                hoxd70_flat[i * 5 + j] = HOXD70[i][j];
        initialized = true;
    }

    return {
        .mat = hoxd70_flat,
        .alphabet_size = 5,
        .gap_open = -GAP_OPEN_PENALTY,
        .gap_extend = -GAP_EXTEND_PENALTY,
        .end_bonus = 0,
        .zdrop = 100,
        .band_width = -1,
        .flag = KSW_EZ_GENERIC_SC
    };
}


Cigar_t globalAlignKSW2(const std::string& ref,
    const std::string& query,
    const KSW2AlignConfig& config) {
    // 编码
    std::vector<uint8_t> ref_enc(ref.size());
    std::vector<uint8_t> query_enc(query.size());
    for (size_t i = 0; i < ref.size(); ++i)
        ref_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(ref[i])];
    for (size_t i = 0; i < query.size(); ++i)
        query_enc[i] = ScoreChar2Idx[static_cast<uint8_t>(query[i])];

    // 比对
    ksw_extz_t ez;
    ksw_reset_extz(&ez);
    ksw_extz2_sse(nullptr,
        static_cast<int>(query_enc.size()), query_enc.data(),
        static_cast<int>(ref_enc.size()), ref_enc.data(),
        config.alphabet_size, config.mat,
        config.gap_open, config.gap_extend,
        config.band_width, config.zdrop, config.end_bonus,
        config.flag, &ez);

    // CIGAR 解析
    Cigar_t cigar;
    for (int i = 0; i < ez.n_cigar; ++i)
        cigar.push_back(ez.cigar[i]);

    if (ez.cigar) free(ez.cigar);
    return cigar;
}


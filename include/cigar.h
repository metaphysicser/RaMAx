#ifndef CIGAR_H
#define CIGAR_H

#include"config.hpp"
#include "alignment/cigar.h"

using Score_t = int_t;

#define MATCH_SCORE_AT  91;
#define MATCH_SCORE_CG  100;
#define MATCH_SCORE_N  -100;

#define GAP_OPEN_PENALTY = -400;
#define GAP_EXTEND_PENALTY = -30;

/// HOXD70 substitution matrix (Chiaromonte et al. 2002)
static const std::array<std::array<int, 5>, 5> HOXD70 = { {
        /*       A      C      G      T     N  */
        /*A*/ {  91,  -114,   -31,  -123, -100 },
        /*C*/ { -114,  100,  -125,   -31, -100 },
        /*G*/ {  -31, -125,   100,  -114, -100 },
        /*T*/ { -123,  -31,  -114,    91, -100 },
        /*N*/ { -100, -100,  -100,  -100, -100}
    } };

static std::array<int8_t, 256> makeCharToIndex() {
    std::array<int8_t, 256> m;
    m.fill(4);
    m['A'] = 0; m['C'] = 1; m['G'] = 2; m['T'] = 3;
    m['a'] = 0; m['c'] = 1; m['g'] = 2; m['t'] = 3;
    return m;
}
static const auto ScoreChar2Idx = makeCharToIndex();

/// 获取 HOXD70 分数
inline int subsScore(char a, char b) {
    return HOXD70[ScoreChar2Idx[a]][ScoreChar2Idx[b]];
}

using CigarUnit = uint32_t;              // 一个 CIGAR 操作
using Cigar_t = std::vector<CigarUnit>; // 整条 CIGAR

// Convert a CIGAR operation and its length to a compact integer representation.
CigarUnit cigarToInt(char operation, uint32_t len);

// Convert a compact integer representation of a CIGAR operation back to its character and length.
void intToCigar(CigarUnit cigar, char& operation, uint32_t& len);

Score_t caculateMatchScore(const char* match, uint_t length);

#endif
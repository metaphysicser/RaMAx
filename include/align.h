#ifndef ALIGN_H
#define ALIGN_H

#include "ksw2.h"
#include "config.hpp"              // 包含基本类型定义，如 int_t、uint_t 等
#include "bindings/cpp/WFAligner.hpp"
extern "C" {
#include "alignment/cigar.h" 
#include "wavefront/wavefront_align.h"
}
// ------------------------------------------------------------------
// 类型定义
// ------------------------------------------------------------------

// 比对得分类型（整数）
using Score_t = int_t;

// ------------------------------------------------------------------
// 基本配对/错配得分设置（宏定义，注意不要加分号）
// ------------------------------------------------------------------
#define MATCH_SCORE_AT  91         // A/T 匹配得分
#define MATCH_SCORE_CG  100        // C/G 匹配得分
#define MATCH_SCORE_N  -100        // 任意与 N 比对惩罚分

#define GAP_OPEN_PENALTY  400     // gap 打开惩罚（较大）
#define GAP_EXTEND_PENALTY  30    // gap 延伸惩罚（较小）


// ------------------------------------------------------------------
// HOXD70 替换矩阵（用于精确比对得分）
// 参考文献：Chiaromonte et al., 2002
// 矩阵顺序：A, C, G, T, N
// 值：匹配为正分，错配为负分，遇到 N 统一惩罚
// ------------------------------------------------------------------
static const std::array<std::array<int, 5>, 5> HOXD70 = { {
        /*       A      C      G      T     N  */
        /*A*/ {  91,  -114,   -31,  -123, -100 },
        /*C*/ { -114,  100,  -125,   -31, -100 },
        /*G*/ {  -31, -125,   100,  -114, -100 },
        /*T*/ { -123,  -31,  -114,    91, -100 },
        /*N*/ { -100, -100,  -100,  -100, -100 }
    } };

// ------------------------------------------------------------------
// 创建字符到 HOXD70 索引的映射
// 支持大小写 ACGT，其他字符默认映射为 N（索引为 4）
// ------------------------------------------------------------------
static std::array<int8_t, 256> makeCharToIndex() {
    std::array<int8_t, 256> m;
    m.fill(4);  // 默认全部映射为 N
    m['A'] = 0; m['C'] = 1; m['G'] = 2; m['T'] = 3;
    m['a'] = 0; m['c'] = 1; m['g'] = 2; m['t'] = 3;
    return m;
}

// 全局常量，用于快速字符查表打分
static const auto ScoreChar2Idx = makeCharToIndex();

// ------------------------------------------------------------------
// 函数：subsScore
// 功能：返回两个碱基字符 a 和 b 在 HOXD70 中的替换得分
// 注意：大小写不敏感，非 ACGT 一律视为 N（惩罚）
// ------------------------------------------------------------------
inline int subsScore(char a, char b) {
    return HOXD70[ScoreChar2Idx[a]][ScoreChar2Idx[b]];
}
// ------------------------------------------------------------------
// 常量：碱基互补表（大写、小写都支持）
// 使用 C++17 constexpr lambda 构建 256 字节映射表
// ------------------------------------------------------------------
inline constexpr std::array<char, 256> BASE_COMPLEMENT = [] {
    std::array<char, 256> m{};
    for (auto& c : m) c = 'N'; // 默认所有字符都映射为 'N'
    m['A'] = 'T';  m['T'] = 'A';
    m['C'] = 'G';  m['G'] = 'C';
    m['a'] = 't';  m['t'] = 'a';
    m['c'] = 'g';  m['g'] = 'c';
    return m;
    }();

inline void reverseComplement(std::string& seq) {
    for (char& c : seq)
        c = BASE_COMPLEMENT[static_cast<unsigned char>(c)];
    std::reverse(seq.begin(), seq.end());
}

// ------------------------------------------------------------------
// CIGAR 表示与转换
// ------------------------------------------------------------------

// 单个 CIGAR 操作的压缩编码（uint32_t）
// 高位表示操作长度，低 4 位为操作类型编码
using CigarUnit = uint32_t;

// 整个 CIGAR 操作序列（压缩形式）
using Cigar_t = std::vector<CigarUnit>;

// ------------------------------------------------------------------
// 函数：cigarToInt
// 功能：将 CIGAR 操作字符（如 'M'）与其长度编码成一个整数
// 编码方式：高 28 位表示长度，低 4 位为操作类型（0=M, 1=I, 等）
// ------------------------------------------------------------------
CigarUnit cigarToInt(char operation, uint32_t len);

// ------------------------------------------------------------------
// 函数：intToCigar
// 功能：将一个压缩整数还原为操作字符与其长度
// 示例：0x50 -> ('M', 5)
// ------------------------------------------------------------------
void intToCigar(CigarUnit cigar, char& operation, uint32_t& len);

// ------------------------------------------------------------------
// 函数：caculateMatchScore
// 功能：对给定的碱基序列（如对齐后片段）按 ACGT 分类打分
// 使用简单模型：A/T 加 MATCH_SCORE_AT，C/G 加 MATCH_SCORE_CG，其他惩罚
// 注：本函数未使用 HOXD70 矩阵
// ------------------------------------------------------------------
Score_t caculateMatchScore(const char* match, uint_t length);

/* ------------------------------------------------------------------
 *  追加/拼接 CIGAR：若两端操作码相同则自动合并
 *  ------------------------------------------------------------------
 *  @param dst  目标 CIGAR（被追加）
 *  @param src  待追加的 CIGAR 片段
 * ------------------------------------------------------------------*/
void appendCigar(Cigar_t& dst, const Cigar_t& src);
/* ------------------------------------------------------------------
 *  追加单个 CIGAR 操作；若与 dst 最后一个操作码一致则合并
 * ------------------------------------------------------------------*/
void appendCigarOp(Cigar_t& dst, char op, uint32_t len);

std::pair<std::string, std::string>
buildAlignment(const std::string& ref_seq,
    const std::string& qry_seq,
    const Cigar_t& cigar);



struct KSW2AlignConfig {
    const int8_t* mat;                  // 一维替换矩阵 (flattened 5x5)
    int alphabet_size;                 // 通常为 5
    int gap_open;                      // gap open penalty (positive)
    int gap_extend;                    // gap extend penalty (positive)
    int end_bonus;                     // 末端奖励分
    int zdrop = 100;                   // Z-drop 剪枝参数
    int band_width = -1;               // -1 表示全矩阵
    int flag = KSW_EZ_GENERIC_SC;      // 默认使用全替换矩阵
};

static int8_t dna5_simd_mat[25];
static void init_simd_mat() {
    static bool done = false;
    if (done) return;
    const int8_t MATCH = 2, MISMATCH = -3, AMBIG = 0;   // N 给 0 分
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            dna5_simd_mat[i * 5 + j] =
            (i == 4 || j == 4) ? AMBIG : (i == j ? MATCH : MISMATCH);
    done = true;
}

//------------------------------------------- 带宽估计
inline int auto_band(int qlen, int tlen,
    double indel_rate = 0.05,
    int    margin = 16)           // 多一点保险
{
    int diff = std::abs(qlen - tlen);
    int extra = static_cast<int>(indel_rate * std::min(qlen, tlen));
    int w = diff + extra + margin;
    return (w + 15) / 16 * 16;                         // 向上取 16 的倍数
}

inline KSW2AlignConfig makeTurboKSW2Config(int qlen, int tlen,
    bool rev_cigar = false,
    double indel_rate = 0.05)
{
    init_simd_mat();                // 确保矩阵已填
    int flags = KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP | KSW_EZ_RIGHT;
    if (rev_cigar) flags |= KSW_EZ_REV_CIGAR;          // 需要时再加
    return {
        .mat = dna5_simd_mat,
        .alphabet_size = 5,
        .gap_open = 8,
        .gap_extend = 1,
        .end_bonus = 0,
        .zdrop = 100,
        .band_width = auto_band(qlen, tlen, indel_rate),
        .flag = flags
    };
}

KSW2AlignConfig makeDefaultKSW2Config();

Cigar_t globalAlignKSW2(const std::string& ref, const std::string& query);
Cigar_t globalAlignWFA2(const std::string& ref, const std::string& query);
#endif



#include "align.h"

// ------------------------------------------------------------------
// 函数：cigarToInt
// 功能：将 CIGAR 操作字符（如 'M'）和其长度编码为一个 32 位无符号整数
// 编码格式：高 28 位表示长度，低 4 位表示操作类型编码（参照 SAM 格式）
// 示例：'M', 5 -> 0x00000050 = (5 << 4) | 0x0
// ------------------------------------------------------------------
uint32_t cigarToInt(char operation, uint32_t len) {
	uint32_t opCode;

	// 将字符型操作转换为对应的数字编码（4 位）
	switch (operation) {
	case 'M': opCode = 0x0; break; // 对齐匹配（可能包含错配）
	case 'I': opCode = 0x1; break; // 插入
	case 'D': opCode = 0x2; break; // 删除
	case '=': opCode = 0x7; break; // 精确匹配（无错配）
	case 'X': opCode = 0x8; break; // 错配
		// TODO: 可以根据需要添加更多操作类型（如 soft clip, hard clip 等）
	default: opCode = 0xF; break; // 未知操作，使用保留码 0xF
	}

	// 将长度左移 4 位，然后用 bitwise OR 合并 opCode
	return (len << 4) | opCode;
}

// ------------------------------------------------------------------
// 函数：intToCigar
// 功能：将压缩后的整数值还原为字符型 CIGAR 操作与其长度
// 参数：
//   - cigar：压缩后的整数（低 4 位为 opCode，其他为长度）
//   - operation：输出的操作字符
//   - len：输出的操作长度
// ------------------------------------------------------------------
void intToCigar(CigarUnit cigar, char& operation, uint32_t& len)
{
	uint32_t opCode = cigar & 0xF;  // 提取低 4 位为操作编码
	len = cigar >> 4;               // 高 28 位表示操作长度

	// 将操作编码转换回对应字符
	switch (opCode) {
	case 0x0: operation = 'M'; break;
	case 0x1: operation = 'I'; break;
	case 0x2: operation = 'D'; break;
	case 0x7: operation = '='; break;
	case 0x8: operation = 'X'; break;
	default:  operation = '?'; break; // 未知操作标记为 '?'
	}
}

// -------------------------------------------------------------------
// 函数：caculateMatchScore
// 功能：计算一段对齐序列的总得分（仅基于字符匹配）
// 参数：
//   - match：表示序列的 char 数组（通常来自 CIGAR 'M' 操作对应的比对段）
//   - length：match 的长度
// 返回：整数得分
// 注：此处未使用 HOXD70 替换矩阵，而是使用简单规则：
//     A/T 匹配 +91，C/G 匹配 +100，其他 -100
// ------------------------------------------------------------------
Score_t caculateMatchScore(const char* match, uint_t length)
{
	Score_t score = 0;

	for (size_t i = 0; i < length; ++i) {
		int id = ScoreChar2Idx[match[i]];  // 获取碱基字符在打分矩阵中的索引

		// 根据碱基种类计算得分
		if (id == 0 || id == 3) {            // A 或 T
			score += MATCH_SCORE_AT;
		}
		else if (id == 1 || id == 2) {       // C 或 G
			score += MATCH_SCORE_CG;
		}
		else {                               // 其他如 N
			score += MATCH_SCORE_N;
		}
	}
	return score;
}

/* ------------------------------------------------------------------
 *  追加/拼接 CIGAR：若两端操作码相同则自动合并
 *  ------------------------------------------------------------------
 *  @param dst  目标 CIGAR（被追加）
 *  @param src  待追加的 CIGAR 片段
 * ------------------------------------------------------------------*/
void appendCigar(Cigar_t& dst, const Cigar_t& src)
{
	if (src.empty()) return;                  // nothing to do
	size_t idx = 0;

	// 1) 检查“拼接点”——dst 的最后一个元素 vs. src 的首元素
	if (!dst.empty()) {
		char op_dst; uint32_t len_dst;
		intToCigar(dst.back(), op_dst, len_dst);

		char op_src; uint32_t len_src;
		intToCigar(src.front(), op_src, len_src);

		// 若操作码相同，则二者合并成一个
		if (op_dst == op_src) {
			dst.back() = cigarToInt(op_dst, len_dst + len_src);
			idx = 1;                         // src 第 0 元素已被合并，后续从 1 开始
		}
	}
	// 2) 将 src 剩余元素逐个 push
	for (; idx < src.size(); ++idx)
		dst.push_back(src[idx]);
}

/* ------------------------------------------------------------------
 *  追加单个 CIGAR 操作；若与 dst 最后一个操作码一致则合并
 * ------------------------------------------------------------------*/
void appendCigarOp(Cigar_t& dst, char op, uint32_t len)
{
	if (len == 0) return;
	if (!dst.empty()) {
		char op_last; uint32_t len_last;
		intToCigar(dst.back(), op_last, len_last);
		if (op_last == op) {                        // 合并
			dst.back() = cigarToInt(op, len_last + len);
			return;
		}
	}
	dst.push_back(cigarToInt(op, len));
}

/* ----------------------------------------------------------
 *  根据 CIGAR 重建 (ref_aln, qry_aln) 两条同宽序列
 *  --------------------------------------------------------*/
std::pair<std::string, std::string>
buildAlignment(const std::string& ref_seq,
	const std::string& qry_seq,
	const Cigar_t& cigar)
{
	std::string ref_aln;  ref_aln.reserve(ref_seq.size() * 1.2);
	std::string qry_aln;  qry_aln.reserve(qry_seq.size() * 1.2);

	size_t i = 0, j = 0;
	for (CigarUnit u : cigar)
	{
		uint32_t len = u >> 4;
		uint8_t  op = u & 0xF;        // 0=M,1=I,2=D,7='=',8='X'

		if (op == 0 || op == 7 || op == 8) {          // match / mismatch
			ref_aln.append(ref_seq.data() + i, len);
			qry_aln.append(qry_seq.data() + j, len);
			i += len;  j += len;
		}
		else if (op == 1) {                       // insertion in query
			ref_aln.append(len, '-');
			qry_aln.append(qry_seq.data() + j, len);
			j += len;
		}
		else if (op == 2) {                       // deletion in query
			ref_aln.append(ref_seq.data() + i, len);
			qry_aln.append(len, '-');
			i += len;
		}
		// 其余 CIGAR 码如 N/H/S 忽略，本例不出现
	}
	return { ref_aln, qry_aln };
}

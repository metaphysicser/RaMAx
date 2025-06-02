/***************************************************************
*  delta2maf.cpp  — convert MUMmer .delta to MAF using actual FASTA
*  Build: g++ -std=c++11 -O2 -o delta2maf delta2maf.cpp
****************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdlib>
#include <algorithm>
#include <iomanip>

/* =============== CLI 参数解析 =============== */
struct Args {
    std::string refFasta;   // 参考 FASTA 文件路径
    std::string qryFasta;   // 查询 FASTA 文件路径
    std::string inDelta;    // 输入 Delta 文件路径
    std::string outMaf;     // 输出 MAF 文件路径
};
void usage(const char* p) {
    std::cerr << "Usage: " << p
              << " --ref /path/to/ref.fa"
              << " --query /path/to/query.fa"
              << " --input file.delta"
              << " --output file.maf\n";
    std::exit(1);
}
Args parse(int ac, char* av[]) {
    Args a;
    for (int i = 1; i < ac; ++i) {
        std::string o(av[i]);
        if      (o == "--ref"    && i + 1 < ac) a.refFasta = av[++i];
        else if (o == "--query"  && i + 1 < ac) a.qryFasta = av[++i];
        else if (o == "--input"  && i + 1 < ac) a.inDelta  = av[++i];
        else if (o == "--output" && i + 1 < ac) a.outMaf   = av[++i];
        else usage(av[0]);
    }
    if (a.refFasta.empty() || a.qryFasta.empty() || a.inDelta.empty() || a.outMaf.empty())
        usage(av[0]);
    return a;
}

/* =============== 反向互补 =============== */
char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default:  return 'N';
    }
}
void revComp(std::string& s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = complement(c);
}

/* =============== 读取 FASTA 到内存 (多染色体支持) =============== */
void loadFasta(const std::string& path, 
               std::vector<std::string>& names, 
               std::vector<std::string>& seqs) 
{
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Cannot open FASTA: " << path << "\n";
        std::exit(1);
    }
    std::string line, curName, curSeq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!curName.empty()) {
                names.push_back(curName);
                seqs.push_back(curSeq);
            }
            // 取 '>' 后第一个空格或换行前为 contig 名
            std::string name = line.substr(1);
            size_t sp = name.find_first_of(" \t");
            if (sp != std::string::npos) name.resize(sp);
            curName = name;
            curSeq.clear();
        } else {
            curSeq += line;
        }
    }
    if (!curName.empty()) {
        names.push_back(curName);
        seqs.push_back(curSeq);
    }
}

/* =============== 最长公共前缀匹配 =============== */
int findLongestPrefixIndex(const std::string& tmpline, 
                           const std::vector<std::string>& names) 
{
    int maxLen = 0, idx = -1;
    for (int i = 0; i < (int)names.size(); ++i) {
        int len = 0;
        int m = std::min((int)tmpline.size(), (int)names[i].size());
        while (len < m && tmpline[len] == names[i][len]) ++len;
        if (len > maxLen) {
            maxLen = len;
            idx = i;
        }
    }
    return idx;
}

/* =============== 读取多行差异列表 =============== */
std::vector<int> parseMultiDiff(std::ifstream& in) {
    std::vector<int> v;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) break;
        std::istringstream ss(line);
        int x;
        while (ss >> x) {
            v.push_back(x);
            if (x == 0) return v;
        }
    }
    return v;
}

void buildAlignment(const std::vector<int>& diffs,
                    const std::string& refSeg,
                    const std::string& qrySeg,
                    std::string& outRef,
                    std::string& outQry)
{
    outRef.clear();  outQry.clear();
    size_t rPos = 0, qPos = 0;

    auto copyMatch = [&](size_t len){
        outRef.append(refSeg, rPos, len);
        outQry.append(qrySeg, qPos, len);
        rPos += len;  qPos += len;
    };

    for (int d : diffs) {
        if (d == 0) break;           // 0 结束一段



        /* 2️⃣ 复制 |d|-1 个匹配（可能为 0） */
        size_t step = std::abs(d) - 1;
        if (step)
            copyMatch(step);

                    /* 1️⃣ 先插入 1 bp 缺口 */
        if (d > 0) {                 // 参考多 1 bp
            outRef.push_back(refSeg[rPos++]);
            outQry.push_back('-');
        } else {                     // 查询多 1 bp
            outRef.push_back('-');
            outQry.push_back(qrySeg[qPos++]);
        }
    }

    /* 3️⃣ 收尾：把剩余片段搬完并给另一端补 '-' */
    if (rPos < refSeg.size()) {
        size_t k = refSeg.size() - rPos;
        outRef.append(refSeg, rPos, k);
    }
    if (qPos < qrySeg.size()) {
        size_t k = qrySeg.size() - qPos;
        outQry.append(qrySeg, qPos, k);
    }
}



/* =============== 主函数 =============== */
int main(int ac, char* av[]) {
    Args args = parse(ac, av);

    // 加载参考与查询 FASTA 到向量
    std::vector<std::string> refNames, refSeqs;
    std::vector<std::string> qryNames, qrySeqs;
    loadFasta(args.refFasta, refNames, refSeqs);
    loadFasta(args.qryFasta, qryNames, qrySeqs);

    std::ifstream in(args.inDelta);
    if (!in) {
        std::cerr << "Cannot open delta file: " << args.inDelta << "\n";
        return 1;
    }
    std::ofstream out(args.outMaf);
    if (!out) {
        std::cerr << "Cannot write MAF file: " << args.outMaf << "\n";
        return 1;
    }

    // 计算列宽，用于对齐输出 contig + 数值
    size_t maxRefName = 0, maxQryName = 0;
    for (auto& nm : refNames) maxRefName = std::max(maxRefName, nm.size());
    for (auto& nm : qryNames) maxQryName = std::max(maxQryName, nm.size());
    size_t colWidth = std::max(maxRefName + 1, maxQryName + 1) + 5;

    // 写 MAF Header
    out << "##maf version=1\n\n";

    std::string line;
    std::string tmpline1, tmpline2;
    int64_t lengthA = 0, lengthB = 0;
    bool inBlock = false;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // 读到 "> refContig qryContig refLen qryLen" 行
            std::istringstream ss(line.substr(1));
            ss >> tmpline1 >> tmpline2 >> lengthA >> lengthB;
            inBlock = true;
            continue;
        }
        if (!inBlock) continue;

        // 读取 coords 行: s1 e1 s2 e2 err idt
        int64_t s1, e1, s2, e2, err, cnt, alnLen;
        std::istringstream ss(line);
        ss >> s1 >> e1 >> s2 >> e2 >> err >> cnt >> alnLen;

        // 读取后续多行差异，直到遇到单行 "0"
        std::vector<int> diffs = parseMultiDiff(in);

        // 在 refNames/qryNames 中寻找最长前缀匹配的索引
        int idxA = findLongestPrefixIndex(tmpline1, refNames);
        int idxB = findLongestPrefixIndex(tmpline2, qryNames);
        if (idxA < 0 || idxB < 0) {
            std::cerr << "ERROR: 无法匹配 contig '" << tmpline1 
                      << "' 或 '" << tmpline2 << "' 至 FASTA 名称\n";
            continue;
        }
        const std::string& refContig = refNames[idxA];
        const std::string& qryContig = qryNames[idxB];
        const std::string& refFull   = refSeqs[idxA];
        const std::string& qryFull   = qrySeqs[idxB];

        // 将 1-based 闭区间转换为 0-based，并检测正／反向
        bool revA = false, revB = false;
        int64_t r1 = s1, r2 = e1;
        int64_t q1 = s2, q2 = e2;
        if (s1 > e1) { std::swap(r1, r2); revA = true; }
        if (s2 > e2) { std::swap(q1, q2); revB = true; }
        int64_t refStart = r1 - 1;
        int64_t refLenSeg = r2 - r1 + 1;
        int64_t qryStart = q1 - 1;
        int64_t qryLenSeg = q2 - q1 + 1;

        // 边界检查
        if (refStart < 0 || refStart + refLenSeg > (int64_t)refFull.size()) {
            std::cerr << "WARNING: ref 片段越界 “" << refContig 
                      << "” : start=" << refStart << ", len=" << refLenSeg 
                      << ", fullLen=" << refFull.size() << "\n";
            continue;
        }
        if (qryStart < 0 || qryStart + qryLenSeg > (int64_t)qryFull.size()) {
            std::cerr << "WARNING: qry 片段越界 “" << qryContig 
                      << "” : start=" << qryStart << ", len=" << qryLenSeg 
                      << ", fullLen=" << qryFull.size() << "\n";
            continue;
        }

        // 提取子串
        std::string refSeg = refFull.substr(refStart, refLenSeg);
        std::string qrySeg = qryFull.substr(qryStart, qryLenSeg);
        if (revA) revComp(refSeg);
        if (revB) revComp(qrySeg);

        // 构建带 gap 的对齐
        std::string gRef, gQry;
        buildAlignment(diffs, refSeg, qrySeg, gRef, gQry);

        // 如果是反向比对，还要修正 MAF 中的显示起点
        int64_t outRefStart = refStart;
        int64_t outQryStart = qryStart;
        char strandA = revA ? '-' : '+';
        char strandB = revB ? '-' : '+';
        if (revA) outRefStart = refFull.size() - (r2);
        if (revB) outQryStart = qryFull.size() - (q2);

        auto countNonGap = [](const std::string& s) {
            return std::count_if(s.begin(), s.end(), [](char c) { return c != '-'; });
        };

        int refNonGapLen = countNonGap(gRef);
        int qryNonGapLen = countNonGap(gQry);

        out << "a score=0\n";
        out << "s " << std::left << std::setw(colWidth) << refContig
            << std::right << std::setw(12) << outRefStart
            << std::setw(12) << refNonGapLen
            << " " << strandA
            << std::setw(12) << (int64_t)refFull.size()
            << " " << gRef << "\n";

        out << "s " << std::left << std::setw(colWidth) << qryContig
            << std::right << std::setw(12) << outQryStart
            << std::setw(12) << qryNonGapLen
            << " " << strandB
            << std::setw(12) << (int64_t)qryFull.size()
            << " " << gQry << "\n\n";

    }

    std::cerr << "MAF written to " << args.outMaf << "\n";
    return 0;
}

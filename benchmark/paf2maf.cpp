// paf2maf.cpp — convert PAF (cg:Z or cs:Z) to MAF using actual FASTA
// Build: g++ -std=c++11 -O2 -o paf2maf paf2maf.cpp
// Author: ChatGPT (adapted from delta2maf style)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cctype>

// =============== CLI 参数解析 ===============
struct Args {
    std::string refFasta;   // reference (target) FASTA
    std::string qryFasta;   // query FASTA
    std::string inPaf;      // input PAF
    std::string outMaf;     // output MAF
};
void usage(const char* p) {
    std::cerr << "Usage: " << p
        << " --ref ref.fa --query query.fa --input file.paf --output out.maf\n";
    std::exit(1);
}
Args parseArgs(int ac, char* av[]) {
    Args a;
    for (int i = 1; i < ac; ++i) {
        std::string o(av[i]);
        if (o == "--ref" && i + 1 < ac) a.refFasta = av[++i];
        else if (o == "--query" && i + 1 < ac) a.qryFasta = av[++i];
        else if (o == "--input" && i + 1 < ac) a.inPaf = av[++i];
        else if (o == "--output" && i + 1 < ac) a.outMaf = av[++i];
        else usage(av[0]);
    }
    if (a.refFasta.empty() || a.qryFasta.empty() || a.inPaf.empty() || a.outMaf.empty())
        usage(av[0]);
    return a;
}

// =============== Reverse-complement ===============
char complement(char c) {
    switch (c) {
    case 'A': return 'T'; case 'C': return 'G';
    case 'G': return 'C'; case 'T': return 'A';
    case 'a': return 't'; case 'c': return 'g';
    case 'g': return 'c'; case 't': return 'a';
    default:  return 'N';
    }
}
void revComp(std::string& s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s) c = complement(c);
}

// =============== FASTA loader ===============
void loadFasta(const std::string& path,
    std::vector<std::string>& names,
    std::vector<std::string>& seqs)
{
    std::ifstream in(path);
    if (!in) { std::cerr << "Cannot open FASTA: " << path << "\n"; std::exit(1); }
    std::string line, curName, curSeq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!curName.empty()) {
                names.push_back(curName);
                seqs.push_back(curSeq);
            }
            curName = line.substr(1);
            auto sp = curName.find_first_of(" \t");
            if (sp != std::string::npos) curName.resize(sp);
            curSeq.clear();
        }
        else {
            curSeq += line;
        }
    }
    if (!curName.empty()) {
        names.push_back(curName);
        seqs.push_back(curSeq);
    }
}

// =============== 找 contig 索引 ===============
int findContig(const std::string& q, const std::vector<std::string>& names) {
    for (size_t i = 0; i < names.size(); ++i) {
        if (q == names[i] || q.find(names[i]) == 0) return (int)i;
    }
    return -1;
}

// =============== Op 结构 ===============
struct Op { char type; int len; }; // 'M'/'I'/'D'

// =============== cg:Z: 解析 ===============
void parseCg(const std::string& cg, std::vector<Op>& ops) {
    std::istringstream ss(cg);
    int v;
    while (ss >> v) {
        if (v == 0) break;
        int m = std::abs(v) - 1;
        if (m > 0) ops.push_back({ 'M', m });
        if (v > 0) ops.push_back({ 'D',1 });
        else     ops.push_back({ 'I',1 });
    }
}

// =============== cs:Z: 解析 ===============
void parseCs(const std::string& cs, std::vector<Op>& ops) {
    int i = 0, n = cs.size();
    while (i < n) {
        char c = cs[i++];
        if (c == '=') {
            int start = i;
            while (i < n && std::isalpha(cs[i])) ++i;
            int L = i - start;
            if (L > 0) ops.push_back({ 'M',L });
        }
        else if (c == '*') {
            // mismatch consumes 1 base each
            if (i + 1 < n && std::isalpha(cs[i]) && std::isalpha(cs[i + 1]))
                i += 2;
            ops.push_back({ 'M',1 });
        }
        else if (c == '+') {
            int start = i;
            while (i < n && std::isalpha(cs[i])) ++i;
            int L = i - start;
            if (L > 0) ops.push_back({ 'I',L });
        }
        else if (c == '-') {
            int start = i;
            while (i < n && std::isalpha(cs[i])) ++i;
            int L = i - start;
            if (L > 0) ops.push_back({ 'D',L });
        }
        else {
            // skip any other
        }
    }
}

// =============== 构建带 gap 的对齐 ===============
void buildAlignment(const std::vector<Op>& ops,
    const std::string& R, const std::string& Q,
    std::string& gR, std::string& gQ)
{
    size_t rpos = 0, qpos = 0;
    gR.clear(); gQ.clear();
    for (auto& o : ops) {
        int L = o.len;
        switch (o.type) {
        case 'M': {
            int a = std::min<size_t>(L, R.size() - rpos);
            a = std::min<size_t>(a, Q.size() - qpos);
            gR.append(R, rpos, a);
            gQ.append(Q, qpos, a);
            rpos += a; qpos += a;
        } break;
        case 'I': {
            int a = std::min<size_t>(L, Q.size() - qpos);
            gR.append(a, '-');
            gQ.append(Q, qpos, a);
            qpos += a;
        } break;
        case 'D': {
            int a = std::min<size_t>(L, R.size() - rpos);
            gR.append(R, rpos, a);
            gQ.append(a, '-');
            rpos += a;
        } break;
        }
    }
    // flush tail
    if (rpos < R.size()) gR.append(R, rpos, R.size() - rpos);
    if (qpos < Q.size()) gQ.append(Q, qpos, Q.size() - qpos);
}

int main(int ac, char* av[]) {
    auto args = parseArgs(ac, av);

    // load FASTA
    std::vector<std::string> refNames, refSeqs, qryNames, qrySeqs;
    loadFasta(args.refFasta, refNames, refSeqs);
    loadFasta(args.qryFasta, qryNames, qrySeqs);

    // open files
    std::ifstream in(args.inPaf);
    std::ofstream out(args.outMaf);
    if (!in) { std::cerr << "Cannot open " << args.inPaf << "\n"; return 1; }
    if (!out) { std::cerr << "Cannot write " << args.outMaf << "\n"; return 1; }

    // col width
    size_t maxR = 0, maxQ = 0;
    for (auto& n : refNames) maxR = std::max(maxR, n.size());
    for (auto& n : qryNames) maxQ = std::max(maxQ, n.size());
    size_t colW = std::max(maxR + 1, maxQ + 1) + 5;

    out << "##maf version=1\n\n";
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);

        // 12 mandatory PAF cols
        std::string qName, tName, strand;
        int64_t qLen, qStart, qEnd;
        int64_t tLen, tStart, tEnd;
        int64_t nMatch, alen; int mapQ;
        ss >> qName >> qLen >> qStart >> qEnd >> strand
            >> tName >> tLen >> tStart >> tEnd
            >> nMatch >> alen >> mapQ;

        // parse tags
        std::string tag, cg, cs;
        while (ss >> tag) {
            if (tag.rfind("cs:Z:", 0) == 0) cs = tag.substr(5);
            else if (tag.rfind("cg:Z:", 0) == 0) cg = tag.substr(5);
        }
        if (cs.empty() && cg.empty()) continue;

        // locate contigs
        int iR = findContig(tName, refNames);
        int iQ = findContig(qName, qryNames);
        bool swapped = false;
        if (iR < 0 || iQ < 0) {
            // 尝试交换：把 Query 当 reference，Target 当 query
            int tmpR = findContig(qName, refNames);
            int tmpQ = findContig(tName, qryNames);
            if (tmpR >= 0 && tmpQ >= 0) {
                iR = tmpR;  iQ = tmpQ;
                std::swap(qName, tName);          // 名称也要交换
                std::swap(qLen, tLen);
                std::swap(qStart, tStart);
                std::swap(qEnd, tEnd);
                strand = (strand == "-") ? "-" : "+";   // strand 在对调后仍保持方向含义
                swapped = true;
            }
            else {
                continue;   // 两种方式都匹配失败，跳过
            }
        }

        auto& Rfull = refSeqs[iR];
        auto& Qfull = qrySeqs[iQ];

        // extract segments
        bool rev = (strand == "-");
        int64_t lenR = tEnd - tStart;
        int64_t lenQ = qEnd - qStart;
        std::string Rseg = Rfull.substr(tStart, lenR);
        std::string Qseg = Qfull.substr(qStart, lenQ);
        if (rev) revComp(Qseg);

        // build ops
        std::vector<Op> ops;
        if (!cs.empty()) parseCs(cs, ops);
        else            parseCg(cg, ops);

        // build gapped sequences
        std::string gR, gQ;
        buildAlignment(ops, Rseg, Qseg, gR, gQ);

        // counts & start pos
        int nonR = std::count_if(gR.begin(), gR.end(), [](char c) {return c != '-';});
        int nonQ = std::count_if(gQ.begin(), gQ.end(), [](char c) {return c != '-';});
        int64_t oR0 = tStart;
        int64_t oQ0 = rev ? (qLen - qEnd) : qStart;

        // emit MAF block
        out << "a score=0\n";
        out << "s " << std::left << std::setw(colW) << tName
            << std::right << std::setw(12) << oR0
            << std::setw(12) << nonR
            << " +" << std::setw(11) << tLen << " " << gR << "\n";
        out << "s " << std::left << std::setw(colW) << qName
            << std::right << std::setw(12) << oQ0
            << std::setw(12) << nonQ
            << " " << (rev ? '-' : '+') << std::setw(12) << qLen << " " << gQ << "\n\n";
    }

    std::cerr << "MAF written to " << args.outMaf << "\n";
    return 0;
}

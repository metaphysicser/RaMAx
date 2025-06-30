// sam2maf.cpp — convert SAM to MAF using actual FASTA
// Build: g++ -std=c++11 -O2 -o sam2maf sam2maf.cpp
// Author: ChatGPT  (2025-06)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cctype>
#include <cstdlib>

/* ==================== CLI ==================== */
struct Args {
    std::string refFasta;   // reference FASTA (chromosomes named in RNAME)
    std::string qryFasta;   // query FASTA (reads/contigs, optional)
    std::string inSam;      // input SAM (plain text)
    std::string outMaf;     // output MAF
};
void usage(const char* p) {
    std::cerr << "Usage: " << p
        << " --ref ref.fa --query query.fa --input aln.sam --output out.maf\n";
    std::exit(1);
}
Args parseArgs(int ac, char* av[]) {
    Args a;
    for (int i = 1;i < ac;++i) {
        std::string o(av[i]);
        if (o == "--ref" && i + 1 < ac) a.refFasta = av[++i];
        else if (o == "--query" && i + 1 < ac) a.qryFasta = av[++i];
        else if (o == "--input" && i + 1 < ac) a.inSam = av[++i];
        else if (o == "--output" && i + 1 < ac) a.outMaf = av[++i];
        else usage(av[0]);
    }
    if (a.refFasta.empty() || a.inSam.empty() || a.outMaf.empty())
        usage(av[0]);
    return a;
}

/* ==================== reverse-complement ==================== */
char comp(char c) {
    switch (c) {
    case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';
    case 'a':return 't';case 'c':return 'g';case 'g':return 'c';case 't':return 'a';
    default:return 'N';
    }
}
void revComp(std::string& s) {
    std::reverse(s.begin(), s.end());
    for (char& c : s)c = comp(c);
}

/* ==================== FASTA loader ==================== */
void loadFasta(const std::string& path,
    std::vector<std::string>& names,
    std::vector<std::string>& seqs) {
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in) { std::cerr << "Cannot open FASTA " << path << "\n"; std::exit(1); }
    std::string line, name, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!name.empty()) { names.push_back(name); seqs.push_back(seq); }
            name = line.substr(1);
            size_t sp = name.find_first_of(" \t");
            if (sp != std::string::npos) name.resize(sp);
            seq.clear();
        }
        else seq += line;
    }
    if (!name.empty()) { names.push_back(name); seqs.push_back(seq); }
}
int findContig(const std::string& n, const std::vector<std::string>& v) {
    for (size_t i = 0;i < v.size();++i)
        if (n == v[i] || n.rfind(v[i], 0) == 0) return (int)i;
    return -1;
}

/* ==================== CIGAR parser ==================== */
struct Op { char type; int len; }; // M / I / D

std::vector<Op> parseCigar(const std::string& cig) {
    std::vector<Op> v;
    long num = 0;
    for (char c : cig) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            num = num * 10 + (c - '0');
        }
        else {
            if (num == 0) { std::cerr << "WARN: zero-len CIGAR op in " << cig << "\n"; }
            switch (c) {
            case 'M': case '=': case 'X': v.push_back({ 'M',(int)num }); break;
            case 'I': v.push_back({ 'I',(int)num }); break;
            case 'D': case 'N': v.push_back({ 'D',(int)num }); break;
            case 'S': case 'H': case 'P': /* clip or pad -> skip */ break;
            default: std::cerr << "WARN: unknown CIGAR op " << c << " ignored\n";
            }
            num = 0;
        }
    }
    return v;
}

/* ==================== alignment builder ==================== */
void buildAlign(const std::vector<Op>& ops,
    const std::string& R, const std::string& Q,
    std::string& gR, std::string& gQ) {
    gR.clear(); gQ.clear();
    size_t r = 0, q = 0;
    for (auto& o : ops) {
        int L = o.len;
        if (o.type == 'M') {
            int a = std::min<int>(L, std::min(R.size() - r, Q.size() - q));
            gR.append(R, r, a); gQ.append(Q, q, a); r += a; q += a;
        }
        else if (o.type == 'I') {
            int a = std::min<int>(L, Q.size() - q);
            gR.append(a, '-'); gQ.append(Q, q, a); q += a;
        }
        else if (o.type == 'D') {
            int a = std::min<int>(L, R.size() - r);
            gR.append(R, r, a); gQ.append(a, '-'); r += a;
        }
    }
    if (r < R.size()) gR.append(R, r, R.size() - r);
    if (q < Q.size()) gQ.append(Q, q, Q.size() - q);
}

/* ==================== main ==================== */
int main(int ac, char* av[]) {
    Args args = parseArgs(ac, av);

    // load FASTA
    std::vector<std::string> refNames, refSeqs, qryNames, qrySeqs;
    loadFasta(args.refFasta, refNames, refSeqs);
    if (!args.qryFasta.empty()) loadFasta(args.qryFasta, qryNames, qrySeqs);

    // open files
    std::ifstream in(args.inSam);
    std::ofstream out(args.outMaf);
    if (!in) { std::cerr << "Cannot open " << args.inSam << "\n"; return 1; }
    if (!out) { std::cerr << "Cannot write " << args.outMaf << "\n"; return 1; }

    // col width
    size_t maxR = 0, maxQ = 0;
    for (auto& n : refNames) maxR = std::max(maxR, n.size());
    for (auto& n : qryNames) maxQ = std::max(maxQ, n.size());
    size_t colW = std::max(maxR + 1, maxQ + 1) + 5;

    out << "##maf version=1\n\n";
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '@') continue;              // header
        std::istringstream ss(line);
        std::string qName, rName, cigar, seq, qual;
        int flag, mapq, pos; std::string dummy;
        ss >> qName >> flag >> rName >> pos >> mapq >> cigar;
        for (int i = 0;i < 2;++i) ss >> dummy;         // RNEXT, PNEXT
        ss >> dummy;                              // TLEN
        ss >> seq >> qual;                          // SEQ QUAL
        if (rName == "*" || cigar == "*") continue;    // unmapped

        // look up reference contig
        int idxR = findContig(rName, refNames);
        if (idxR < 0) { std::cerr << "WARN: reference contig " << rName << " not found\n"; continue; }
        const std::string& Rfull = refSeqs[idxR];

        // query sequence
        std::string Qfull;
        if (!args.qryFasta.empty()) {
            int idxQ = findContig(qName, qryNames);
            if (idxQ >= 0) Qfull = qrySeqs[idxQ];
            else Qfull = seq;                 // fallback to inline SEQ
        }
        else Qfull = seq;

        int readLen = Qfull.size();
        bool rev = (flag & 16);

        // compute reference slice
        int64_t refStart0 = pos - 1;            // SAM is 1-based
        std::vector<Op> ops = parseCigar(cigar);

        // derive aligned lengths (non-clip)
        int refAligned = 0, qryAligned = 0;
        for (auto& o : ops) {
            if (o.type == 'M') { refAligned += o.len; qryAligned += o.len; }
            else if (o.type == 'I') qryAligned += o.len;
            else if (o.type == 'D') refAligned += o.len;
        }
        std::string Rseg = Rfull.substr(refStart0, refAligned);

        // derive query slice (remove leading soft clips)
        int leadClip = 0, trailClip = 0;
        long num = 0; bool leading = true;
        for (char c : cigar) {
            if (std::isdigit(static_cast<unsigned char>(c))) { num = num * 10 + (c - '0'); }
            else {
                if (c == 'S') { if (leading) leadClip += num; else trailClip += num; }
                else leading = false;
                num = 0;
            }
        }
        std::string Qseg = Qfull.substr(leadClip, qryAligned);
        if (rev) revComp(Qseg);

        // build gapped alignment
        std::string gR, gQ;
        buildAlign(ops, Rseg, Qseg, gR, gQ);

        // MAF start
        int qryStart = rev ? (readLen - trailClip - qryAligned) : leadClip;
        int nonR = std::count_if(gR.begin(), gR.end(), [](char c) {return c != '-';});
        int nonQ = std::count_if(gQ.begin(), gQ.end(), [](char c) {return c != '-';});

        // output
        out << "a score=" << mapq << "\n";
        out << "s " << std::left << std::setw(colW) << rName
            << std::right << std::setw(12) << refStart0
            << std::setw(12) << nonR << " +"
            << std::setw(11) << Rfull.size() << " " << gR << "\n";

        out << "s " << std::left << std::setw(colW) << qName
            << std::right << std::setw(12) << qryStart
            << std::setw(12) << nonQ << " "
            << (rev ? '-' : '+') << std::setw(12) << readLen << " " << gQ << "\n\n";
    }

    std::cerr << "MAF written to " << args.outMaf << "\n";
    return 0;
}

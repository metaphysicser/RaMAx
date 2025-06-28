// paf2maf.cpp — convert PAF to MAF using actual FASTA (reference = target, query = query)
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

/* ================= Command‑line ================= */
struct Args {
    std::string refFasta;   // target FASTA (reference)
    std::string qryFasta;   // query FASTA
    std::string inPaf;      // input PAF
    std::string outMaf;     // output MAF
};

void usage(const char* p) {
    std::cerr << "Usage: " << p
        << " --ref ref.fa --query query.fa --input aln.paf --output out.maf\n";
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

/* ================= Reverse‑complement ================= */
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

/* ================= FASTA loader ================= */
void loadFasta(const std::string& path,
    std::vector<std::string>& names,
    std::vector<std::string>& seqs) {
    std::ifstream in(path);
    if (!in) { std::cerr << "Cannot open FASTA: " << path << "\n"; std::exit(1); }
    std::string line, curName, curSeq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!curName.empty()) { names.push_back(curName); seqs.push_back(curSeq); }
            curName = line.substr(1);
            size_t sp = curName.find_first_of(" \t");
            if (sp != std::string::npos) curName.resize(sp);
            curSeq.clear();
        }
        else curSeq += line;
    }
    if (!curName.empty()) { names.push_back(curName); seqs.push_back(curSeq); }
}

/* ================= CIGAR parser ================= */
using Op = std::pair<char, int>; // (op, length)

std::vector<Op> parseCigar(const std::string& cigar) {
    std::vector<Op> v;
    int num = 0;
    for (char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            num = num * 10 + (c - '0');
        }
        else {
            v.emplace_back(c, num);
            num = 0;
        }
    }
    return v;
}

void buildAlignment(const std::vector<Op>& ops,
    const std::string& refSeg,
    const std::string& qrySeg,
    std::string& gRef,
    std::string& gQry) {
    size_t r = 0, q = 0;
    gRef.clear(); gQry.clear();
    for (auto& kv : ops) {
        char op = kv.first;
        int  len = kv.second;
        switch (op) {
        case '=': case 'X':
            gRef.append(refSeg.substr(r, len));
            gQry.append(qrySeg.substr(q, len));
            r += len; q += len;
            break;
        case 'I': // insertion to ref (bases in query only)
            gRef.append(len, '-');
            gQry.append(qrySeg.substr(q, len));
            q += len;
            break;
        case 'D': // deletion from ref (bases in ref only)
            gRef.append(refSeg.substr(r, len));
            gQry.append(len, '-');
            r += len;
            break;
        }
    }
}

/* ================= Main ================= */
int main(int ac, char* av[]) {
    Args args = parseArgs(ac, av);

    // Load FASTA sequences
    std::vector<std::string> refNames, refSeqs;
    std::vector<std::string> qryNames, qrySeqs;
    loadFasta(args.refFasta, refNames, refSeqs);
    loadFasta(args.qryFasta, qryNames, qrySeqs);

    // Prepare output
    std::ifstream in(args.inPaf);
    if (!in) { std::cerr << "Cannot open PAF: " << args.inPaf << "\n"; return 1; }
    std::ofstream out(args.outMaf);
    if (!out) { std::cerr << "Cannot write MAF: " << args.outMaf << "\n"; return 1; }

    // Column width for pretty MAF
    size_t maxRef = 0, maxQry = 0;
    for (auto& n : refNames) maxRef = std::max(maxRef, n.size());
    for (auto& n : qryNames) maxQry = std::max(maxQry, n.size());
    size_t colW = std::max(maxRef + 1, maxQry + 1) + 5;

    // MAF header
    out << "##maf version=1\n\n";

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);

        // PAF mandatory fields (query first, target second)
        std::string qName, tName, strand; // strand is query relative to target
        int64_t qLen, qStart, qEnd;
        int64_t tLen, tStart, tEnd;
        int64_t nMatch, alnLen; int mapQ;

        ss >> qName >> qLen >> qStart >> qEnd >> strand
            >> tName >> tLen >> tStart >> tEnd
            >> nMatch >> alnLen >> mapQ; // consume MAPQ too

        if (!ss) continue; // malformed

        // Grab cg:Z: tag (CIGAR) if present
        std::string tag, cg;
        while (ss >> tag) {
            if (tag.rfind("cg:Z:", 0) == 0) {
                cg = tag.substr(5); // after cg:Z:
                break;
            }
        }
        if (cg.empty()) continue; // no detailed CIGAR -> skip

        // Locate sequences in FASTA vectors
        auto itRef = std::find_if(refNames.begin(), refNames.end(), [&](const std::string& n) {return n == tName || tName.rfind(n, 0) == 0;});
        auto itQry = std::find_if(qryNames.begin(), qryNames.end(), [&](const std::string& n) {return n == qName || qName.rfind(n, 0) == 0;});
        if (itRef == refNames.end() || itQry == qryNames.end()) {
            std::cerr << "WARNING: contig not found for line: " << qName << " -> " << tName << "\n";
            continue;
        }
        size_t idxRef = std::distance(refNames.begin(), itRef);
        size_t idxQry = std::distance(qryNames.begin(), itQry);
        const std::string& refFull = refSeqs[idxRef];
        const std::string& qryFull = qrySeqs[idxQry];

        // Compute segment boundaries (0‑based, half‑open)
        bool revQry = (strand == "-");
        int64_t refSegLen = tEnd - tStart;
        int64_t qrySegLen = qEnd - qStart;

        if (refSegLen <= 0 || qrySegLen <= 0) continue;
        if (tStart <0 || tEnd > tLen || qStart<0 || qEnd>qLen) continue;
        if (tEnd > (int64_t)refFull.size() || qEnd > (int64_t)qryFull.size()) {
            std::cerr << "WARNING: coordinates exceed FASTA length at " << qName << " / " << tName << "\n";
            continue;
        }

        std::string refSeg = refFull.substr(tStart, refSegLen);
        std::string qrySeg = qryFull.substr(qStart, qrySegLen);
        if (revQry) revComp(qrySeg);

        // Build gapped alignment strings
        auto ops = parseCigar(cg);
        std::string gRef, gQry;
        buildAlignment(ops, refSeg, qrySeg, gRef, gQry);

        int refNonGap = std::count_if(gRef.begin(), gRef.end(), [](char c) {return c != '-';});
        int qryNonGap = std::count_if(gQry.begin(), gQry.end(), [](char c) {return c != '-';});

        // Compute MAF start positions (0‑based)
        int64_t outRefStart = tStart;
        int64_t outQryStart = revQry ? (qLen - qEnd) : qStart;

        // Write block
        out << "a score=0\n";
        out << "s " << std::left << std::setw(colW) << tName
            << std::right << std::setw(12) << outRefStart
            << std::setw(12) << refNonGap
            << " +" << std::setw(11) << tLen  // plus one space already in +
            << " " << gRef << "\n";

        out << "s " << std::left << std::setw(colW) << qName
            << std::right << std::setw(12) << outQryStart
            << std::setw(12) << qryNonGap
            << " " << (revQry ? '-' : '+')
            << std::setw(12) << qLen
            << " " << gQry << "\n\n";
    }

    std::cerr << "MAF written to " << args.outMaf << "\n";
    return 0;
}

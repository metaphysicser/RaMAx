// filter_maf_by_species.cpp
// ------------------------------------------------------------------
// Build: g++ -std=c++11 -O3 -fopenmp -o filter_maf_by_species filter_maf_by_species.cpp
// ------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <cstdlib>
#include <omp.h>

/* =============== CLI ================= */
struct Args {
    std::string in, out, ref, query;
    int threads = 1;
};
void usage(const char* p) {
    std::cerr << "Usage: " << p
              << " -i <input.maf> -o <output.maf> "
              << "--ref <species> --query <species> "
              << "[-t <threads>]\n";
    std::exit(EXIT_FAILURE);
}
Args parseArgs(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string opt(argv[i]);
        if ((opt=="-i"||opt=="--input") && i+1<argc)       a.in = argv[++i];
        else if ((opt=="-o"||opt=="--output") && i+1<argc) a.out = argv[++i];
        else if (opt=="--ref"   && i+1<argc)               a.ref = argv[++i];
        else if (opt=="--query" && i+1<argc)               a.query = argv[++i];
        else if ((opt=="-t"||opt=="--threads") && i+1<argc)a.threads = std::atoi(argv[++i]);
        else usage(argv[0]);
    }
    if(a.in.empty()||a.out.empty()||a.ref.empty()||a.query.empty()) usage(argv[0]);
    if(a.threads<=0) a.threads=1;
    return a;
}

/* =============== helpers ================= */
using Block = std::vector<std::string>;

inline std::string speciesName(const std::string& sLine) {
    // s <src> ...
    std::istringstream ss(sLine);
    std::string tag, src;
    ss >> tag >> src;
    return src.substr(0, src.find('.'));
}

struct SLine {
    std::string tag, src, startStr, sizeStr, strand, srcSize, seq;
};

SLine parseSLine(const std::string& line) {
    std::istringstream ss(line);
    SLine s;
    ss >> s.tag >> s.src >> s.startStr >> s.sizeStr >> s.strand >> s.srcSize;
    std::getline(ss, s.seq);           // leading space included
    if (!s.seq.empty() && s.seq[0]==' ') s.seq.erase(0,1);   // trim one leading space
    return s;
}

std::string buildSLine(const SLine& s) {
    std::ostringstream os;
    os << s.tag << ' ' << s.src << ' '
       << s.startStr << ' ' << s.sizeStr << ' '
       << s.strand   << ' ' << s.srcSize << ' '
       << s.seq << '\n';
    return os.str();
}

/* =============== read blocks + capture header ================= */
std::vector<Block> readBlocks(const std::string& path, std::string& preamble) {
    std::ifstream in(path.c_str());
    if(!in){ std::cerr<<"Cannot open "<<path<<"\n"; std::exit(1);}

    std::vector<Block> blocks;
    Block blk;
    std::string line;
    bool seenFirstA = false;

    while(std::getline(in,line)) {
        line.push_back('\n');

        if(!seenFirstA && line.rfind("a ",0)!=0) {      // 还没到第一个 block
            preamble += line;                           // 记录 header 行
            continue;
        }

        if(line.rfind("a ",0)==0){                      // new header
            seenFirstA = true;
            if(!blk.empty()) blocks.emplace_back(std::move(blk));
            blk.clear();
            blk.emplace_back(std::move(line));
        } else if(line=="\n") {                         // blank line -> end block
            if(!blk.empty()) blocks.emplace_back(std::move(blk));
            blk.clear();
        } else {
            blk.emplace_back(std::move(line));
        }
    }
    if(!blk.empty()) blocks.emplace_back(std::move(blk));
    return blocks;
}


/* =============== core block processing ================= */
std::string processBlock(Block& blk,
                         const std::unordered_set<std::string>& keep)
{
    const std::string& header = blk.front();
    std::vector<size_t> chosenIdx;   // index of s-lines kept
    for(size_t i=1;i<blk.size();++i){
        if(blk[i].size() && blk[i][0]=='s'){
            if(keep.count(speciesName(blk[i])))
                chosenIdx.push_back(i);
        }
    }
    if(chosenIdx.size()!=2) return "";              // 必须同时含两物种

    /* 解析两条 s 行 */
    SLine A = parseSLine(blk[chosenIdx[0]]);
    SLine B = parseSLine(blk[chosenIdx[1]]);
    if(A.seq.size()!=B.seq.size()) return "";       // 块有问题

    /* 删除双 gap 列 */
    std::string newA, newB;
    newA.reserve(A.seq.size()); newB.reserve(B.seq.size());
    for(size_t i=0;i<A.seq.size();++i){
        if(A.seq[i]=='-' && B.seq[i]=='-') continue; // 双 gap 列删除
        newA.push_back(A.seq[i]);
        newB.push_back(B.seq[i]);
    }
    if(newA.empty()) return "";                     // 删除后空

    /* 重新计算 size 字段 */
    auto ungap = [](const std::string& s){
        size_t n=0; for(char c:s) if(c!='-') ++n; return n;
    };
    A.seq = newA;  A.sizeStr = std::to_string(ungap(newA));
    B.seq = newB;  B.sizeStr = std::to_string(ungap(newB));

    /* 输出新的 block */
    std::ostringstream os;
    os << header
       << buildSLine(A)
       << buildSLine(B)
       << '\n';
    return os.str();
}

/* =============== main ================= */
int main(int argc,char* argv[]){
    Args args = parseArgs(argc, argv);
    omp_set_num_threads(args.threads);

    std::string preamble;                                     // <- NEW
    std::vector<Block> blocks = readBlocks(args.in, preamble);
    std::vector<std::string> outputs(blocks.size());

    std::unordered_set<std::string> keep{args.ref, args.query};

    #pragma omp parallel for schedule(dynamic)
    for(size_t i=0;i<blocks.size();++i){
        outputs[i] = processBlock(blocks[i], keep);
    }

    std::ofstream out(args.out.c_str());
    if(!out){ std::cerr<<"Cannot write "<<args.out<<"\n"; std::exit(1);}

    out << preamble;                                          // <- NEW
    for(const auto& s : outputs) if(!s.empty()) out<<s;

    std::cerr << "Filtered MAF saved to: " << args.out << "\n";
    return 0;
}

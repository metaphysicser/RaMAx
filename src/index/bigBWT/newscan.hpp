#ifndef NEWSCAN_HPP
#define NEWSCAN_HPP

extern "C" {
#include "xerrors.h"
}
#include <vector>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
}
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <pthread.h>
// KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace __gnu_cxx;
extern pthread_mutex_t map_mutex;

#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// ========================== Word Stats and Command Args ===================
struct word_stats {
    string str;
    occ_int_t occ;
    word_int_t rank = 0;
};

struct Args {
    string inputFileName = "";
    size_t w = 10;            // sliding window size and its default
    size_t p = 100;           // modulus for establishing stopping w-tuples
    bool SAinfo = false;   // compute SA information
    bool is_fasta = false;  // read a fasta file
    bool compress = false; // parsing called in compress mode
    int th = 0;              // number of helper threads
    int verbose = 0;         // verbosity level
};

// ========================== KR Window Structure ============================
struct KR_window {
    int wsize;
    int* window;
    int asize;
    const uint64_t prime = 1999999973;
    uint64_t hash;
    uint64_t tot_char;
    uint64_t asize_pot;

    KR_window(int w) : wsize(w) {
        asize = 256;
        asize_pot = 1;
        for (int i = 1;i < wsize;i++)
            asize_pot = (asize_pot * asize) % prime;
        window = new int[wsize];
        reset();
    }

    void reset() {
        for (int i = 0;i < wsize;i++) window[i] = 0;
        hash = tot_char = 0;
    }

    uint64_t addchar(int c) {
        int k = tot_char++ % wsize;
        hash += (prime - (window[k] * asize_pot) % prime);
        hash = (asize * hash + c) % prime;
        window[k] = c;
        return hash;
    }

    string get_window() {
        string w = "";
        int k = (tot_char - 1) % wsize;
        for (int i = k + 1;i < k + 1 + wsize;i++)
            w.append(1, window[i % wsize]);
        return w;
    }

    ~KR_window() {
        delete[] window;
    }
};

uint64_t kr_hash(string s);
void save_update_word(string& w, unsigned int minsize, map<uint64_t, word_stats>& freq, FILE* tmp_parse_file, FILE* last, FILE* sa, uint64_t& pos);
uint64_t process_file(Args& arg, map<uint64_t, word_stats>& wordFreq);
void writeDictOcc(Args& arg, map<uint64_t, word_stats>& wfreq, vector<const string*>& sortedDict);
void remapParse(Args& arg, map<uint64_t, word_stats>& wfreq);
bool is_gzipped(std::string fname);

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end, start, end; // input
  size_t skipped, parsed, words;  // output
  FILE *parse, *last, *sa;
} mt_data;


void* mt_parse(void* dx);
uint64_t mt_process_file(Args& arg, std::map<uint64_t, word_stats>& wf);
void* mt_parse_fasta(void* dx);
uint64_t mt_process_file_fasta(Args& arg, std::map<uint64_t, word_stats>& wf);
#endif
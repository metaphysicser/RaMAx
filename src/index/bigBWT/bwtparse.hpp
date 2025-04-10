#ifndef BWT_PARSE_HPP
#define BWT_PARSE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include "gsa/gsacak.h" // Assuming this is the correct path to this file
#include "utils.h"       // Assuming this is the correct path to this file

namespace BWTParse {

    // Typedef for the SA index (used for both 32-bit and 64-bit)
    typedef uint_t sa_index_t;

    // Command-line arguments structure
    struct Args {
        char* filename;  // Base name for the files
        bool SAinfo;           // Flag to indicate if SA info is required
        int th;                // Number of segments for last and sa files
    };

    // Read the parse file and add a 0 EOS symbol at the end
    static uint32_t* read_parse(char* filename, long* tsize)
    {
        FILE* parse = open_aux_file(filename, EXTPARSE, "rb");
        // get file size
        if (fseek(parse, 0, SEEK_END) != 0) die("parse fseek");
        long nn = ftell(parse);
        // check input file is OK
        if (nn % 4 != 0) {
            printf("Invalid input file: size not multiple of 4\n");
            exit(1);
        }
#if !M64
        // if in 32 bit mode, the number of words is at most 2^31-2
        if (nn / 4 > 0x7FFFFFFE) {
            printf("Input containing more than 2^31-2 phrases!\n");
            printf("Please use 64 bit version\n");
            exit(1);
        }
#else
        // if in 64 bit mode, the number of words is at most 2^32-2 (for now)
        if (nn / 4 > 0xFFFFFFFEu) {
            printf("Input containing more than 2^32-2 phrases!\n");
            printf("This is currently a hard limit\n");
            exit(1);
        }
#endif
        // printf("Parse file contains %ld words\n", nn / 4);
        long n = nn / 4;
        // ------ allocate and read text file, len is n+1 for the EOS 
        uint32_t* Text = (uint32_t*)malloc((n + 1) * sizeof(*Text));  // 强制类型转换

        if (Text == NULL) die("malloc failed (Text)");
        rewind(parse);

        // read the array in one shot
        assert(sizeof(*Text) == 4);
        size_t s = fread(Text, sizeof(*Text), n, parse);
        if (s != n) {
            char* msg = NULL;
            int e = asprintf(&msg, "read parse error: %zu vs %ld\n", s, n);
            (void)e; die(msg);
        }
        if (fclose(parse) != 0) die("parse file close");
        Text[n] = 0; // sacak needs a 0 eos 
        *tsize = n;
        return Text;
    }


    // Function to compute the SA using the sacak algorithm
    //static sa_index_t* compute_SA(uint32_t* Text, long n, long k) {
    //    sa_index_t* SA = new sa_index_t[n];
    //    int depth = sacak_int(Text, SA, n, k);
    //    if (depth >= 0)
    //        std::cout << "SA computed with depth: " << depth << std::endl;
    //    else
    //        throw std::runtime_error("Error computing SA");

    //    return SA;
    //}

    // Function to load the .last file into an array
    static uint8_t* load_last(Args* arg, long n)
    {
        // open .last file for reading
        mFile* lastin = mopen_aux_file(arg->filename, EXTLST, arg->th);
        // allocate and load the last array
        uint8_t* last = (uint8_t*)malloc(n);
        if (last == NULL) die("malloc failed (LAST)");
        size_t s = mfread(last, 1, n, lastin);
        if (s != n) die("last read");
        if (mfclose(lastin) != 0) die("last file close");
        return last;
    }

    // Function to load the .sa info file into an array (if required)
    static uint8_t* load_sa_info(const Args* arg, long n) {
        if (!arg->SAinfo) return nullptr;

        mFile* fin = mopen_aux_file(arg->filename, EXTSAI, arg->th);
        uint8_t* sai = new uint8_t[n * IBYTES];
        size_t s = mfread(sai, IBYTES, n, fin);
        if (s != n)
            throw std::runtime_error("Error reading SA info file");
        mfclose(fin);
        return sai;
    }

    // Function to open the output file for SA information
    static FILE* open_sa_out(const Args* arg) {
        if (!arg->SAinfo) return nullptr;
        return open_aux_file(arg->filename, EXTBWSAI, "wb");
    }


} // namespace BWTParse

#endif // BWT_PARSE_HPP

//
// Created by Matt Vincent on 5/19/15.
//

#ifndef KALLISTO_SHOWALIGNMENTS_H
#define KALLISTO_SHOWALIGNMENTS_H

#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "MinCollector.h"

#include "common.h"

void printStringVector(const std::vector<std::string>& v, std::ostream& o) {
    int i = 0;
    for (auto x : v) {
        if (i>0) {
            o << ", ";
        }
        o << x;
        i++;
    }
}

template<typename Index, typename TranscriptCollector>
void ShowAlignments(Index& index, const ProgramOptions& opt, TranscriptCollector& tc) {
    // need to receive an index map
    std::ios_base::sync_with_stdio(false);

    int tlencount = 10000;
    size_t numreads = 0;
    size_t nummapped = 0;

    bool paired = !opt.single_end;

    gzFile fp1 = 0, fp2 = 0;
    kseq_t *seq1 = 0, *seq2 = 0;
    std::vector<std::pair<int,int>> v1, v2;
    v1.reserve(1000);
    v2.reserve(1000);

    int l1 = 0,l2 = 0; // length of read

    if (paired) {
        std::cerr << "[quant] running in paired-end mode" << std::endl;
    } else {
        std::cerr << "[quant] running in single-end mode" << std::endl;
    }

    for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {
        if (paired) {
            std::cerr << "[quant] will process pair " << (i/2 +1) << ": "  << opt.files[i] << std::endl
            << "                             " << opt.files[i+1] << std::endl;
        } else {
            std::cerr << "[quant] will process file " << i+1 << ": " << opt.files[i] << std::endl;
        }
    }

    // for each file
    std::cerr << "[quant] finding pseudoalignments for the reads ..." << std::endl;
    std::cerr.flush();

    for (int i = 0; i < opt.files.size(); i += (paired) ? 2 : 1) {

        fp1 = gzopen(opt.files[i].c_str(), "r");
        seq1 = kseq_init(fp1);
        if (paired) {
            fp2 = gzopen(opt.files[i+1].c_str(),"r");
            seq2 = kseq_init(fp2);
        }



        // for each read
        while (true) {
            l1 = kseq_read(seq1);
            if (paired) {
                l2 = kseq_read(seq2);
            }
            if (l1 <= 0) {
                break;
            }
            if (paired && l2 <= 0) {
                break;
            }

            numreads++;
            v1.clear();
            v2.clear();
            // process read
            index.match(seq1->seq.s, seq1->seq.l, v1);
            if (paired) {
                index.match(seq2->seq.s, seq2->seq.l, v2);
            }

            // collect the target information
            int ec = tc.collect(v1, v2, !paired);
            if (ec != -1) {
                nummapped++;
                std::cout << seq1->seq.s << "\t";
                const std::vector<int> &vec = index.ecmap[ec];
                std::vector<std::string> target_names;
                for (int a : vec) {
                    target_names.push_back(index.target_names_[a]);
                }
                printStringVector(target_names, std::cout);
                std::cout << std::endl;
            } else {
                std::cout << seq1->seq.s << "\t" << std::endl;
            }

            if (paired && 0 <= ec &&  ec < index.num_trans && tlencount > 0) {
                //bool allSame = true;
                bool allSame = (v1[0].first == ec && v2[0].first == ec) && (v1[0].second == 0 && v2[0].second == 0);

                if (allSame) {
                    // try to map the reads
                    int tl = index.mapPair(seq1->seq.s, seq1->seq.l, seq2->seq.s, seq2->seq.l, ec);
                    if (0 < tl && tl < tc.flens.size()) {
                        tc.flens[tl]++;
                        tlencount--;
                    }
                }
            }

        }
        gzclose(fp1);
        if (paired) {
            gzclose(fp2);
        }
    }

    kseq_destroy(seq1);
    if (paired) {
        kseq_destroy(seq2);
    }

    std::cerr << " done" << std::endl;

    //std::cout << "betterCount = " << betterCount << ", out of betterCand = " << betterCand << std::endl;

    std::cerr << "[quant] processed " << numreads << " reads, " << nummapped << " reads pseudoaligned" << std::endl;

    // write output to outdir
    if (opt.write_index) {
        std::string outfile = opt.output + "/counts.txt";
        std::ofstream of;
        of.open(outfile.c_str(), std::ios::out);
        tc.write(of);
        of.close();
    }
}


#endif //KALLISTO_SHOWALIGNMENTS_H

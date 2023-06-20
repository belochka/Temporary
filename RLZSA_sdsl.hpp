#ifndef _RLZSA_sdsl_H_
#define _RLZSA_sdsl_H_

#include <stdio.h>
#include <stdlib.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <chrono>
#include <ctime>
#include <ratio>
#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

using namespace std;
using namespace std::chrono;

class RLZSA_sdsl {

public:
    RLZSA_sdsl() {
    }

    ~RLZSA_sdsl() {
    }

    void load(char *rlzfilename, char *safilename) {
        int junk;

        //read in the original SA
        ifstream safin(safilename, ios::in | ios::binary);
        safin.seekg(0, safin.end);
        _n = safin.tellg() / sizeof(uint32_t);
        safin.seekg(0);
        _sa = new uint32_t[_n];
        safin.read((char *) _sa, _n * sizeof(uint32_t));
        safin.close();

        ifstream fin(rlzfilename, ios::in | ios::binary);
        fin.read((char *) &_r, sizeof(size_t));
        uint32_t * reference = new uint32_t[_r];
        fin.read((char *) reference, _r * sizeof(uint32_t));
        _reference = sdsl::int_vector<>(_r, 0, sdsl::bits::hi(_n) + 2);
        std::copy(reference, reference + _r, _reference.begin());
        delete[] reference;

        fin.read((char *) &_z, sizeof(size_t));
        uint32_t * phrases = new uint32_t[_z];
        fin.read((char *) phrases, _z * sizeof(uint32_t));
        _phrasesIntWidth = sdsl::bits::hi(_n - 1) + 1;
        _phrases = sdsl::int_vector<>(_z, 0, _phrasesIntWidth + 1);
        std::copy(phrases, phrases + _z, _phrases.begin());
        delete[] phrases;

        fin.read((char *) &_z, sizeof(size_t));
        uint32_t *starts = new uint32_t[_z + 1];
        fin.read((char *) starts, _z * sizeof(uint32_t));
        //put a dummy on the end of starts
        starts[_z] = _n;
        fin.close();
        predecessorDS = sdsl::sd_vector<>(starts, starts + _z);

        //loop through and replace literal _phrases with their original SA values
        for (uint32_t i = 0; i < _z; i++) {
            if (_phrases[i] == starts[i]) {
                //literal phrase
                _phrases[i] = _sa[starts[i]] | (1 << _phrasesIntWidth);
            }
        }

////        cerr << "About to compute partial sums and phrase lengths\n";
//        //build predecessor data structure for phrase starting positions
////        _numpblocks = (_z / _pbs);
////        if (_z / _pbs) {
////            _numpblocks++;
////        }
////        _partialSums = new uint32_t[_numpblocks + 1];
////        _partialSums[_numpblocks] = _n;
//        _phraseLengths = new uint16_t[_z];
//        for (int i = 0; i < _z; i++) {
////            if (i % _pbs == 0) {
////                _partialSums[i / _pbs] = starts[i];
////            }
//            _phraseLengths[i] = starts[i + 1] - starts[i];
//        }
    }

    void queryOneCell(char *filename) {
        //read in pattern ranges
        ifstream fin(filename, ios::in | ios::binary);
        //fin.open(rangeFilename, ios::in | ios::binary);
        fin.seekg(0, fin.end);
        size_t ncells = fin.tellg() / sizeof(uint32_t);
        fin.seekg(0);
        vector <uint32_t> cells;
        for (int i = 0; i < ncells; i++) {
            uint32_t position;
            fin.read((char *) &position, sizeof(uint32_t));
//         cerr << position << '\n';
            cells.emplace_back(position);
        }
        fin.close();

        uint64_t search;
#ifdef DCC2000
        uint64_t sumSearchTimes = 0, timesToRun = 10;
        for (int runi = 0; runi < timesToRun; runi++) {
#endif
        auto t1 = high_resolution_clock::now();
        for (int i = 0; i < ncells; i++) {
#ifdef perocc
            auto tOcc1 = high_resolution_clock::now();
#endif
            uint32_t occ;
            uint32_t offset = 0;
            uint32_t p;
            uint32_t sum = biggestLowerOrEqual(predecessorDS, cells[i], p);
            const uint32_t phraseLength = getPhraseLength(predecessorDS, sum, p);
            if ((_phrases[p] & (1 << _phrasesIntWidth))) {
                // if it is a literal phrase
                occ = _phrases[p] & ((1 << _phrasesIntWidth) - 1);
            } else {
                //repeat phrase
//                if (p == 0) {
//                    cerr << "ALARM! The very first phrase isn't a literal phrase, cannot proceed!\n";
//                    exit(1);
//                }
                uint32_t prevSaValue = _phrases[p - 1] & ((1 << _phrasesIntWidth) - 1);
                offset = cells[i] - sum;    // adjustment for starting position inside first phrase
                for (uint32_t j = 0; j < offset; j++) {
                    prevSaValue = _reference[_phrases[p] + j] + prevSaValue - _n;
                }
                occ = _reference[_phrases[p] + offset] + prevSaValue - _n;
            }
#ifdef perocc
            auto tOcc2 = high_resolution_clock::now();
            uint64_t cellTime = std::chrono::duration_cast<std::chrono::nanoseconds>(tOcc2 - tOcc1).count();
            fprintf(stdout, "%d,%d,%u,%lu\n", runi, i, offset, cellTime);
#endif
            if (occ != _sa[cells[i]]) {
                cout << "AHTUNG!!! occ != _sa[" << cells[i] << "]\n" << occ << " instead of " << _sa[cells[i]] << '\n';
                exit(1);
            }
        }

        auto t2 = high_resolution_clock::now();
        search = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
#ifdef DCC2000
        sumSearchTimes += search;
        }
        search = (uint64_t)((double)sumSearchTimes / timesToRun);
#endif
    }

    void query(char *rangeFilename) {
        //read in pattern ranges
        ifstream fin(rangeFilename, ios::in | ios::binary);
        //fin.open(rangeFilename, ios::in | ios::binary);
        fin.seekg(0, fin.end);
        size_t nranges = fin.tellg() / sizeof(uint32_t) / 2;
        fin.seekg(0);
        vector <pair<uint32_t, uint32_t>> ranges;
        for (int i = 0; i < nranges; i++) {
            uint32_t first, second;
            fin.read((char *) &first, sizeof(uint32_t));
            fin.read((char *) &second, sizeof(uint32_t));
//         cerr << first << ' ' << second << '\n';
            ranges.emplace_back(first, second);
        }
        fin.close();

//        cerr << "Ranges are in memory\n";
//        cerr << "_n = " << _n << "\n";
//        cerr << "_reference.size() = " << _reference.size() << "\n";
//        cerr << "_phrases.size() = " << _phrases.size() << "\n";
//        cerr << "_phraseLengths.size() = " << _z << "\n";

        uint64_t totalNumberOfOccurrences = 0;
//        uint64_t totalThatIsSupposedToBeHere = 0;

        //extract occurrences from RLZSA_sdsl
        uint64_t search;
#ifdef DCC2000
        uint64_t sumSearchTimes = 0, timesToRun = 10;
    for (int runi = 0; runi < timesToRun; runi++) {
#endif
        auto t1 = high_resolution_clock::now();
        for (int i = 0; i < nranges; i++) {
#ifdef perocc
            auto tOcc1 = high_resolution_clock::now();
#endif
            uint32_t nocc = ranges[i].second - ranges[i].first + 1;
            uint32_t *occs = new uint32_t[nocc];
//            totalThatIsSupposedToBeHere += nocc;
//            //uint32_t firstPhrase = pred(ranges[i].first);
////            cerr << "nocc: " << nocc << '\n';
////            cerr << "ranges[" << i << "].first: " << ranges[i].first << " ranges[" << i << "].second: " << ranges[i].second << '\n';
            uint32_t p;
            uint32_t sum = biggestLowerOrEqual(predecessorDS, ranges[i].first, p);
            uint32_t o = 0;
            //deal with first phrase
            const uint32_t phraseLength = getPhraseLength(predecessorDS, sum, p);
//            cerr << "sum: " << sum << " p: " << p << " len: " << phraseLength << " o: " << o << " _phrasesIntWidth: " << _phrasesIntWidth << " _n: " << _n << '\n';
            if ((_phrases[p] & (1 << _phrasesIntWidth))) {
                // if it is a literal phrase
                occs[o++] = _phrases[p] & ((1 << _phrasesIntWidth) - 1);
            } else {
                //repeat phrase
//                if (p == 0) {
//                    cerr << "ALARM! The very first phrase isn't a literal phrase, cannot proceed!\n";
//                    exit(1);
//                }
                uint32_t prevSaValue = _phrases[p - 1] & ((1 << _phrasesIntWidth) - 1);
                uint32_t offset = ranges[i].first - sum;    // adjustment for starting position inside first phrase
                uint32_t end = phraseLength - offset;
                if (end > nocc) { end = nocc; }
                for (uint32_t j = 0; j < offset; j++) {
                    prevSaValue = _reference[_phrases[p] + j] + prevSaValue - _n;
                }
                for (uint32_t j = offset; j < offset + end; j++) {
                    occs[o++] = _reference[_phrases[p] + j] + prevSaValue - _n;
                    prevSaValue = occs[o - 1];
                }
            }
//            cerr << "o (after first phrase): " << o << '\n';
            sum += phraseLength;
            p++;
            //deal with "middle" _phrases
            while (sum <= ranges[i].second) {
                const uint32_t phraseLength = getPhraseLength(predecessorDS, sum, p);
//                cerr << "sum: " << sum << " p: " << p << " len: " << phraseLength << " o: " << o << '\n';
                uint32_t len = phraseLength;
                if (sum + len > ranges[i].second) {
                    //we are in the last phrase, adjust len accordingly
                    len = ranges[i].second - sum + 1;
                }
                if ((_phrases[p] & (1 << _phrasesIntWidth))) {
                    //literal phrase
                    occs[o++] = _phrases[p] & ((1 << _phrasesIntWidth) - 1);
                } else {
                    //repeat phrase
//                    if (p == 0) {
//                        cerr << "ALARM! The very first phrase (p == 0) isn't a literal phrase, cannot proceed!\n";
//                        exit(1);
//                    }
                    uint32_t prevSaValue = _phrases[p - 1] & ((1 << _phrasesIntWidth) - 1);
                    for (uint32_t j = 0; j < len; j++) {
                        occs[o++] = _reference[_phrases[p] + j] + prevSaValue - _n;
                        prevSaValue = occs[o - 1];
                    }
                }
                p++;
                sum += len;
            }
#ifdef perocc
            auto tOcc2 = high_resolution_clock::now();
            uint64_t patternTime = std::chrono::duration_cast<std::chrono::nanoseconds>(tOcc2 - tOcc1).count();
            fprintf(stdout, "%d,%d,%u,%lu\n", runi, i, nocc, patternTime);
#endif

//            cerr << "sum: " << sum << " p: " << p << '\n';
//            cerr << "o (at end): " << o << '\n';
            totalNumberOfOccurrences += o;
//                        for (uint32_t oi = 0; oi < o; oi++) {
//                            if (occs[oi] != _sa[ranges[i].first + oi]) {
//                                fprintf(stderr, "AHTUNG! occs[%u] = %u, should be SA[%u] = %u\n", oi, occs[oi], ranges[i].first + oi, _sa[ranges[i].first + oi]);
//                                exit(1);
//                            }
//                        }
            delete[] occs;
        }

        auto t2 = high_resolution_clock::now();
        search = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
#ifdef DCC2000
        sumSearchTimes += search;
    }
    search = (uint64_t)((double)sumSearchTimes / timesToRun);
#endif
        //cerr << "number of patterns n = " << n << endl;
        //cerr << "pattern length m = " << m << endl;
        //cerr << "total number of occurrences  occ_t = " << occ_tot << endl;

//        cout << "Total time : " << search << " milliseconds" << endl;
//        cout << "Total number of occurrences : " << totalNumberOfOccurrences << endl;
//        if (totalNumberOfOccurrences != totalThatIsSupposedToBeHere) {
//            cout << "AHTUNG!!! totalNumberOfOccurrences != totalThatIsSupposedToBeHere\n" << totalNumberOfOccurrences << " instead of " << totalThatIsSupposedToBeHere << '\n';
//            exit(1);
//        }

    	cout << search << "," << totalNumberOfOccurrences << ",";

    }

    //using std::chrono::high_resolution_clock;
    //using std::chrono::duration_cast;
    //using std::chrono::duration;


    //uint64_t *locate(uint64_t sp, uint64_t ep){
    //   if(sp>ep) return nullptr;
    //   uint64_t occ = ep-sp+1;
    //   uint64_t results = new uint64_t[occ];
    //   uint64_t offset = 0; //how far into the first phrase we start
    //   uint64_t phraseIndex = predecessor(sp,offset); //index into array of <pos,abs> tuples; sets offset
    //   uint64_t phraseLength = getLength(phraseIndex);
    //   do{
    //      uint64_t pos = _posabs[phraseIndex].first;
    //      uint64_t abs = _posabs[phraseIndex].second;
    //      for(uint64_t j=offset; j<phraseLength-offset; j++){
    //         results[ri++] = abs + _reference[pos+j]; //could replace with a memcpy (and SIMD?), maybe
    //      }
    //      phraseIndex++;
    //      offset = 0;
    //   }while(ri < occ);
    //   return results;
    //}

private:

    uint32_t biggestLowerOrEqual(const sdsl::sd_vector<> & sdVector, const uint32_t v, uint32_t &endingPosLow) {
        size_t groupOfTheAnswer;

        uint32_t group = (v >> sdVector.wl);
        uint32_t groupBitHighPos = sdVector.high_0_select(group + 1);
        endingPosLow = groupBitHighPos - group - 1;

        bool found = false;
        if (sdVector.high[groupBitHighPos - 1] == 1) {
            uint32_t prevGroupStartingPosHigh = group > 0 ? sdVector.high_0_select(group) : 0;
            size_t startingPosLow = prevGroupStartingPosHigh - (group - 1);
            if ((group << sdVector.wl) + sdVector.low[startingPosLow] <= v) {
                for (size_t i = startingPosLow + 1; i <= endingPosLow; i++) {
                    if ((group << sdVector.wl) + sdVector.low[i] > v) {
                        endingPosLow = i - 1;
                        groupOfTheAnswer = group;
                        found = true;
                        break;
                    }
                }
            } else {
                endingPosLow = prevGroupStartingPosHigh - (group - 1) - 1;
            }
        }

        if (!found) {
            groupOfTheAnswer = sdVector.high_1_select(endingPosLow + 1) + 1 - endingPosLow - 1;
        }

        return (groupOfTheAnswer << sdVector.wl) + sdVector.low[endingPosLow];
    }

    uint32_t getPhraseLength(const sdsl::sd_vector<> & sdVector, const uint32_t v, uint32_t endingPosLow) {
        uint32_t group = (v >> sdVector.wl);
        if (sdVector.high[group + endingPosLow + 1] == 1) {
            return sdVector.low[endingPosLow + 1] - sdVector.low[endingPosLow];
        }

        uint32_t newGroup;
//        if (sdVector.high[group + endingPosLow + 2] == 1) {
//            newGroup = group + 1;
//        } else {
            uint64_t sdVecSize = sdVector.size();
            uint64_t next1StartingPosHigh;

            //  find startingPositionLow of the current group using 64-bit word checks
            next1StartingPosHigh = group + endingPosLow + 2;
            uint64_t len = (next1StartingPosHigh + 64 >= sdVecSize) ? sdVecSize - next1StartingPosHigh : 64;
            uint64_t word;
            while (next1StartingPosHigh < sdVecSize
                   && ((word = sdVector.high.get_int(next1StartingPosHigh, len)) == ((uint64_t) 0))) {
                next1StartingPosHigh += len;
                len = (next1StartingPosHigh + 64 >= sdVecSize) ? sdVecSize - next1StartingPosHigh : 64;
            }
            uint64_t index = (uint64_t) __builtin_ctzl((word << (64 - len)));
            next1StartingPosHigh += index;

//            uint32_t next1StartingPosHigh = sdVector.high_1_select(endingPosLow + 2);
            newGroup = next1StartingPosHigh - endingPosLow - 1;
//        }

        return ((newGroup << sdVector.wl) + sdVector.low[endingPosLow + 1]) - ((group << sdVector.wl) + sdVector.low[endingPosLow]);
    }

    size_t _n;   //length of SA and text
    size_t _z;   //number of _phrases
    size_t _r;   //number of symbols in _reference
    size_t _phrasesIntWidth;

//    uint16_t *_phraseLengths;
//    uint32_t _pbs = 64; //predecessor block size
//    uint32_t _numpblocks = 64; //predecessor block size
//    uint32_t *_partialSums; //every _pbs^th
    sdsl::int_vector<> _reference;
    sdsl::int_vector<> _phrases;
    uint32_t *_sa;

    sdsl::sd_vector<> predecessorDS;
};

#endif

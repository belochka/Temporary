
#ifndef AXFXLX_NEB_SEB_HPP
#define AXFXLX_NEB_SEB_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <cmath>        // for exp2()
#include "select_support_intv.hpp"
#include "x86intrin.h"

using namespace sdsl;

template<
        class t_rank_1     = typename bit_vector::rank_1_type>
class rank_support_nebseb;

template<
        class t_rank_1     = typename bit_vector::rank_1_type>
class Neb_Seb {
public:
    typedef bit_vector::size_type                   size_type;

private:
    size_type m_size = 0;           // length of the original bit vector
    uint8_t   m_wl   = 0;           // log n - log m, where n is the length of the original bit vector

    int_vector<>          m_low;    // vector for the least significant bits of the positions of the m ones
    bit_vector            m_neb;    // vector that says which buckets are empty and which are not (NEB - non empty buckets)
//    bit_vector            m_seb;  // vector that shows in unary the size of each non-empty bucket (SEB - size of existing bucket). Code a bucket of size k in unary as "0^(k-1)1"
    int_vector<>          m_seb;    // vector that shows in unary the size of each non-empty bucket (SEB - size of existing bucket). In decimal format

    t_rank_1              m_neb_1_rank;          // rank support for the ones in m_neb
    select_support_intv   m_seb_1_select;

public:
    const uint8_t&               wl             = m_wl;
    const bit_vector&            neb            = m_neb;
    const int_vector<>&          seb            = m_seb;    // NB! remember that values stores = (value - 1), and add 1!!
    const int_vector<>&          low            = m_low;

    uint32_t _partialBlockSize = 16;  //predecessor block size
    uint32_t _logPartialBlockSize = 4;  //predecessor block size
    size_type *_partialSums;                 //every _partialBlockSize^th. last one = sum of all
    size_type _numblocks;

    const t_rank_1&               neb_1_rank    = m_neb_1_rank;
    const select_support_intv&    seb_select    = m_seb_1_select; // NB! You better not touch that if m_wl > 8! I mean it.

public:
    Neb_Seb() {
    }

    Neb_Seb(const Neb_Seb& nebseb)
    {
        copy(nebseb);
    }

    Neb_Seb(Neb_Seb&& nebseb)
    {
        *this = std::move(nebseb);
    }

    Neb_Seb(const bit_vector& bv)
    {
        m_size = bv.size();
        size_type _n = util::cnt_one_bits(bv);
        uint8_t logn = bits::hi(_n)+1;
        uint8_t logm = bits::hi(m_size)+1;
        if (logn == logm) {
            --logn;    // to ensure logn-logm > 0
        }
        m_wl    = logm - logn;
        m_low = int_vector<>(_n, 0, m_wl);
        size_type maxBucketSize = exp2(m_wl);
        m_neb = bit_vector(ceil(m_size / maxBucketSize) + 1);
        m_seb = int_vector<>(_n, 0UL, m_wl);

//        std::cout << "m_wl = " << std::to_string(m_wl) << std::endl;
//        std::cout << "64 / m_wl = " << std::to_string(64 / m_wl) << std::endl;

        const size_type* bvp = bv.data();
        size_type sebi = 0, bucketSize = 0;
        for (size_type i = 0, lowi = 0, last_high = 0; i < (bv.size() + 63) / 64; ++i, ++bvp) {
            size_type position = 64*i;
            size_type  w = *bvp;
            while (w) {  // process bit_vector word by word
                uint8_t offset = bits::lo(w);
                w >>= offset;   // note:  w >>= (offset+1) can not be applied for offset=63!
                position += offset;
                if (position >= bv.size()) // check that we have not reached the end of the bitvector
                    break;

                // (1) handle high part
                size_type cur_high = position >> m_wl;
                if (cur_high != last_high && bucketSize != 0) {
                    m_seb[sebi++] = bucketSize - 1;
                    bucketSize = 0;
                }
                last_high = cur_high;

                // (2) handle low part
                m_neb[cur_high] = 1;          // write 1 for the entry
                m_low[lowi++] = position;       // "int_vector truncates the most significant logm bits" (c) Simon Gog in sd_vector.hpp
                bucketSize += 1;
                position += 1;
                w >>= 1;
            }
        }
        if (bucketSize != 0) {
            m_seb[sebi++] = bucketSize - 1;
        }
        m_seb.resize(sebi);

        util::init_support(m_neb_1_rank, &m_neb);

        if (m_wl <= 8) {
            util::init_support(m_seb_1_select, &m_seb);
        } else {
            // compute _partialSums
//            _partialBlockSize = 64 / m_wl;
            _numblocks = ((m_seb.size() - 1) >> _logPartialBlockSize) + 1;
            _partialSums = new size_type[_numblocks];
            size_type sum = 0;
            size_type mask = (1 << _logPartialBlockSize) - 1; // '1111'
            for (size_type i = 0; i < m_seb.size(); i++) {
                if (i != 0 && (i % _partialBlockSize == 0)) {
                    if ( (i & mask) != 0 ) {
                        fprintf(stderr, "Error 1: i = %lu, (i & mask) = %lu, _logPartialBlockSize = %lu\n", i, (size_type)(i & mask), (size_type)_logPartialBlockSize);
                        exit(1);
                    }
                }
                if ((i != 0) && ((i & mask) == 0)) {
                    _partialSums[(i >> _logPartialBlockSize) - 1] = sum;
                }
                sum += m_seb[i] + 1;
            }
            _partialSums[_numblocks - 1] = _n;
        }

//        std::cout << "Construction is done" << std::endl;
    }

    //! Returns the size of the original bit vector.
    size_type size()const
    {
        return m_size;
    }

    ~Neb_Seb() = default;

    Neb_Seb& operator=(const Neb_Seb& v)
    {
        if (this != &v) {
            copy(v);
        }
        return *this;
    }

    Neb_Seb& operator=(Neb_Seb&& v)
    {
        if (this != &v) {
            m_size = v.m_size;
            m_wl   = v.m_wl;
            m_low  = std::move(v.m_low);
            m_neb  = std::move(v.m_neb);
            m_seb  = std::move(v.m_seb);

            m_neb_1_rank = std::move(v.m_neb_1_rank);
            m_neb_1_rank.set_vector(&m_neb);

            if (m_wl <= 8) {
                util::init_support(m_seb_1_select, &m_seb);
//                m_seb_1_select = std::move(v.m_seb_1_select);  // .width() is 64 here => assertion fails
//                m_seb_1_select.set_vector(&m_seb);
            } else {
                _partialBlockSize = v._partialBlockSize;
                _partialSums = v._partialSums;
                _numblocks = v._numblocks;
            }
        }
        return *this;
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "size");
        written_bytes += write_member(m_wl, out, child, "wl");
        written_bytes += m_low.serialize(out, child, "low");
        written_bytes += m_neb.serialize(out, child, "neb");
        written_bytes += m_seb.serialize(out, child, "seb");

        // not necessary to be serialized, can be calculated from previous
        written_bytes += m_neb_1_rank.serialize(out, child, "neb_1_rank");
        if (m_wl <= 8) {
            written_bytes += m_seb_1_select.serialize(out, child, "seb_1_select");
        } else {
            written_bytes += write_member(_partialBlockSize, out, child, "partialBlockSize");
            written_bytes += write_member(_numblocks, out, child, "numblocks");
            written_bytes += _numblocks * sizeof(_partialSums[0]);      // Don't want to write serialization of an array
        }

        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream& in)
    {
        read_member(m_size, in);
        read_member(m_wl, in);
        m_low.load(in);
        m_neb.load(in);
        m_seb.load(in);

        // can be calculated from previous
        m_neb_1_rank.load(in, &m_neb);
        if (m_wl <= 8) {
            m_seb_1_select.load(in, &m_seb);
        } else {
            read_member(_partialBlockSize, in);
            read_member(_numblocks, in);

//        // compute _partialSums
//        _partialBlockSize = 64 / m_wl;
//        _numblocks = ((m_seb.size() - 1) / _partialBlockSize) + 1;
            _partialSums = new size_type[_numblocks];
            size_type sum = 0;
            size_type mask = (1 << _logPartialBlockSize) - 1; // '1111'
            for (size_type i = 0; i < m_seb.size(); i++) {
                if ((i != 0) && ((i & mask) == 0)) {
                    _partialSums[(i >> _logPartialBlockSize) - 1] = sum;
                }
                sum += m_seb[i] + 1;
            }
            _partialSums[_numblocks - 1] = sum;
        }
    }


private:
    void copy(const Neb_Seb& v)
    {
        m_size = v.m_size;
        m_wl   = v.m_wl;
        m_low  = v.m_low;
        m_neb = v.m_neb;
        m_seb = v.m_seb;

        m_neb_1_rank = v.m_neb_1_rank;
        m_neb_1_rank.set_vector(&m_neb);

        if (m_wl <= 8) {
            m_seb_1_select = v.m_seb_1_select;
            m_seb_1_select.set_vector(&m_seb);
        } else {
            _partialBlockSize = v._partialBlockSize;
            _partialSums = v._partialSums;
            _numblocks = v._numblocks;
        }
    }
};


//! Rank data structure for NEBSEB + sd_vector.low
template<class t_rank_1>
class rank_support_nebseb {
public:
    typedef bit_vector::size_type size_type;
    typedef Neb_Seb<t_rank_1> bit_vector_type;
private:
    const bit_vector_type* m_v;

    uint8_t m_bucketSearchType;    // Type = 0 (linear) if m_wl <= 3, type = 1 (broadword) if m_wl in [4,8], and type = 2 (binary) if m_wl > 8.
    size_type m_pdep_mask;

public:

    explicit rank_support_nebseb(const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type rank(size_type i)const
    {
        assert(m_v != nullptr);
        assert(i <= m_v->size());

        size_type result = 0;
        size_type bucket = i >> m_v->wl;
        size_type notEmptyBuckets = m_v->neb_1_rank(bucket);

        if (notEmptyBuckets >= m_v->seb.size()) {
            result = (m_v->wl > 8) ? m_v->_partialSums[m_v->_numblocks - 1] : m_v->seb_select(notEmptyBuckets);
        } else {
            size_type startFrom = 0;
            if (m_v->wl <= 8) {
                result = m_v->seb_select(notEmptyBuckets);
            } else {
                if (notEmptyBuckets > m_v->_partialBlockSize) {
                    size_type i = notEmptyBuckets >> m_v->_logPartialBlockSize;
                    result = m_v->_partialSums[i - 1];
                    startFrom = i << m_v->_logPartialBlockSize;
                }
                for (size_type i = startFrom; i < notEmptyBuckets; i++) {
                    result += m_v->seb[i] + 1;
                }
            }

            if (m_v->neb[bucket]) { // if NEB[.] == 1 => bucket not empty => search in bucket as well. Bucket in low in positions [result, result + SEB[notEmptyBuckets]-1]
                size_type val_low = i & bits::lo_set[m_v->wl];
                startFrom = result;

                switch (m_bucketSearchType) {
                    case 0:
                        search_linear(startFrom, startFrom + m_v->seb[notEmptyBuckets], val_low, result);
                        break;
                    case 1:
                        search_broadword(startFrom, startFrom + m_v->seb[notEmptyBuckets], val_low, result);
                        break;
                    case 2:
                        search_binary(startFrom, startFrom + m_v->seb[notEmptyBuckets], val_low, result);
                        break;
                }
            }
        }

        return result;
    }

    size_type operator()(size_type i)const
    {
        return rank(i);
    }

    size_type size()const
    {
        return m_v->size();
    }

    void set_vector(const bit_vector_type* v=nullptr)
    {
        m_v = v;
        initData();
    }

    rank_support_nebseb& operator=(const rank_support_nebseb& rs)
    {
        if (this != &rs) {
            set_vector(rs.m_v);
            m_bucketSearchType = rs.m_bucketSearchType;
        }
        return *this;
    }

    void swap(rank_support_nebseb& v) {
        if (this != &v) {
            std::swap(m_bucketSearchType, v.m_bucketSearchType);
        }
    }

    void load(std::istream&, const bit_vector_type* v=nullptr)
    {
        set_vector(v);
    }

    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
    {
        return serialize_empty_object(out, v, name, this);
    }

private:
    void initData() {
        if (m_v != nullptr) {
            m_pdep_mask = ((1 << m_v->wl) - 1) * (~0ULL/255);

//            // Type = 0 (linear) if m_wl <= 3, type = 1 (broadword) if m_wl in [4,8], and type = 2 (binary) if m_wl > 8.
            if (m_v->wl <= 3) {
                m_bucketSearchType = 0;
            } else if (m_v->wl >= 4 && m_v->wl <= 8) {
                m_bucketSearchType = 1;
            } else {
                m_bucketSearchType = 2;
            }

//            if (m_v->wl <= 4) {
//                m_bucketSearchType = 0;
//            } else {
//                m_bucketSearchType = 2;
//            }
        }
    }

    void copy(const rank_support_nebseb& v)
    {
        m_bucketSearchType = v.m_bucketSearchType;
    }




    inline size_type search_linear(size_type startFrom, size_type endIncluding, size_type val_low, size_type& result)const
    {
        for (size_type i = startFrom; i <= endIncluding; i++) {
            if (m_v->low[i] < val_low) {
                result++;
            } else {
                break;
            }
        }
    }

////                    // (3) Rajeev's binsrch
    inline size_type search_binary(size_type startFrom, size_type endIncluding, size_type val_low, size_type& result)const
    {
        if (m_v->low[endIncluding] >= val_low) {
            size_type s = startFrom, e = endIncluding;
            size_type mid;
            while (e - s > 0) {
                mid = (s + e) >> 1;
                if (m_v->low[mid] < val_low)
                    s = mid + 1;
                else
                    e = mid;

            }
            result += s - startFrom;
        } else {
            result += endIncluding - startFrom + 1;
        }
    }

    inline size_type search_broadword(size_type startFrom, size_type endIncluding, size_type val_low, size_type& result)const
    {
        /* This test is needed -- not just an optimization. See below. */
        if (val_low == 0) return 0;

        size_t bitstoget = 8 * m_v->wl;
        __m64 query = _mm_set1_pi8(val_low);

        while (startFrom + 7 <= endIncluding) {
            result += avx_rank(m_v->low.get_int(startFrom * m_v->wl, bitstoget), query, m_pdep_mask);
            startFrom += 8;
        }

        /* the word that is read using get_int below only contains
      (e - s + 1) real values, not 8.  The other values will be
      taken as 0 by avx_rank and will artificially increase the
      value returned by avx_rank when val_low > 0.  This needs
      to be corrected. */

        size_t correction = 8 - (endIncluding - startFrom + 1);

        if(startFrom <= endIncluding) {
            uint8_t remaining = (endIncluding - startFrom + 1) * m_v->wl;
            result += avx_rank(m_v->low.get_int(startFrom * m_v->wl, remaining), query, m_pdep_mask) - correction;
        }

        return 0;
    }

    inline int avx_rank(const uint64_t x, const __m64 y, uint64_t pdep_mask) const {
        uint64_t e_x = _pdep_u64(x, pdep_mask);
        __m64 __me_x = *(__m64*) &e_x;

        /* Find all positions in __me_x that are >= key */
        __m64 cmpge = _m_pcmpeqb (__me_x, _mm_max_pu8 (y, __me_x));

        return (8 - (_popcnt64((int64_t) cmpge)>>3));
    }
};

#endif //AXFXLX_NEB_SEB_HPP

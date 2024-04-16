/*
 PepinInt

 Copyright (c) 2021, Mate Soos and Kuldeep S. Meel. All rights reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#pragma once

#include <vector>
#include <cstdint>
#include <cassert>
#include <limits>
#include <random>
#include <iostream>
#include <gmp.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include "pepin.h"

using std::vector;
using std::cout;
using std::endl;


#define release_assert(a) \
    do { \
        if (!(a)) {\
            fprintf(stderr, "*** ASSERTION FAILURE in %s() [%s:%d]: %s\n", \
            __FUNCTION__, __FILE__, __LINE__, #a); \
            abort(); \
        } \
    } while (0)

using PepinNS::Lit;
namespace PepinIntNS {

typedef uint16_t WeightPrec;


#define print_verb(n, X) \
    if (verbosity >= n) cout << "c [dnfs] " << X << endl

struct Weight {
    Weight() {}
    Weight(WeightPrec _dividend, WeightPrec _divisor) :
        dividend(_dividend),
        divisor(_divisor)
    {
        assert(dividend <= divisor);
        bits_of_entropy = calc_bits_of_entropy();
    }
    WeightPrec dividend = 1;
    WeightPrec divisor = 2;
    double bits_of_entropy = 1.0;

    double calc_bits_of_entropy() const {
        if (dividend == 1 && divisor == 2) return 1;
        // TODO must calculate a sane entropy value. This is a VERY
        // weak estimate now
        return 0;
    }
};

typedef unsigned char value;

struct Sample {
    uint32_t sol_at;
    uint32_t ws_at;
    uint32_t gen_point;

    bool equals_sol_w(const Sample& other, uint32_t nvars,
                      const vector<value>& sol, const vector<WeightPrec>& w) const {
        for(uint32_t i = 0; i < nvars; i ++) {
            if (sol[sol_at+i] != sol[other.sol_at+i]) {
                return false;
            }
        }
        for(uint32_t i = 0; i < nvars; i ++) {
            if (w[ws_at+i] != w[other.ws_at+i]) {
                return false;
            }
        }
        return true;
    }
};

class Elems {
public:
    ~Elems() {
        free(data);
    }

    uint64_t size() const {
        return num_elems;
    }

    value operator[](const uint64_t at) const {
        uint64_t at_pos = at/4; //4 values per byte
        uint32_t sub_at = (at % 4)<<1;
        return (data[at_pos] >> sub_at) & 3;
    }

    void set(uint64_t at, value val) {
        assert(val <= 3);
        uint64_t at_pos = at/4; //4 values per bit
        uint32_t sub_at = (at % 4)*2;

        //Clear value
        data[at_pos] &= ~((uint8_t)3 << sub_at);

        //Set value
        val <<= sub_at;
        data[at_pos] |= val;
    }

    void fill_unset(uint64_t at, uint64_t num) {
        assert(num % 4 == 0);
        assert(at % 4 == 0);
        uint64_t at_pos = at/4; //4* values per bit
        memset(data+at_pos, 0xff, num/4);
    }

    void insert_unset(uint64_t num) {
        assert(num % 4 == 0);

        //we align it to pagesize
        const auto pagesize = getpagesize();
        uint8_t* data2;
        auto ret = posix_memalign((void**)&data2, pagesize, curr_size+num/4);
        assert(ret == 0);
        if (data) memcpy(data2, data, curr_size);
        memset(data2+curr_size, 0xff, num/4);
        curr_size+=num/4;
        free(data);
        data = data2;

        //and tell the kernel that we'll use it by accessing it randomly
        madvise(data, curr_size, MADV_RANDOM);
        num_elems += num;
    }


private:
    uint64_t curr_size = 0;
    uint8_t* data = NULL;
    uint64_t num_elems = 0;

};

class ElemDat {
public:
    explicit ElemDat (uint32_t  _dnf_cl_num)
        : dnf_cl_num(_dnf_cl_num)
    {}

    bool empty = false;
    uint32_t checked = 0;
    uint32_t updated = 0;
    uint64_t dnf_cl_num = 0; //step at which we generated it
};

class Bucket {
public:
    Bucket (const bool& _all_default_weights, const uint32_t& _verbosity) :
        all_default_weights(_all_default_weights),
        verbosity(_verbosity)
    {}

    uint64_t add_lazy_common(const uint64_t dnf_cl_num) {
        uint64_t at;
        if (empties.empty()) {
            elems.insert_unset(nvars);
            at = elems.size()-nvars;
            elems_dat.push_back(ElemDat(dnf_cl_num));
        } else {
            const size_t orig_at = empties.back();
            at = orig_at*nvars;
            empties.pop_back();
            elems.fill_unset(at, nvars);
            elems_dat[orig_at] = ElemDat(dnf_cl_num);
        }
        return at;
    }

    void add(const value* sol, const uint64_t dnf_cl_num) {
        assert(nvars > 0);
        assert((elems.size() % nvars) == 0);
        assert(elems_dat.size() == elems.size()/nvars);

        const uint64_t at = add_lazy_common(dnf_cl_num);
        for(uint32_t i = 0; i < nvars; i++) elems.set(at+i, sol[i]);
        size++;
    }

    void add_lazy(const vector<Lit>& cl, const uint64_t dnf_cl_num);

    uint64_t get_size() const {
        return size;
    }

    size_t nVars() const {
        return nvars;
    }

    void remove_sol(const vector<Lit>& cl,
                    const vector<Weight>& var_weights,
                    std::mt19937_64& mtrand);

    void remove(const size_t at) {
        assert(elems_dat[at].empty == false);
        assert(elems.size()/nvars > at);

        empties.push_back(at);
        elems_dat[at].empty = true;
        size--;
    }

    void remove_half(std::mt19937_64& mtrand);
    void print_contents() const;
    void set_nvars(uint64_t _nvars) {
        assert(_nvars > 0);
        assert(nvars == 0);
        nvars = _nvars;
    }

    void print_elems_stats(const uint64_t tot_num_dnf_cls) const;

private:
    uint64_t size = 0;
    uint32_t nvars = 0;
    //vector<value> elems;
    Elems elems;
    vector<uint32_t> empties;
    vector<ElemDat> elems_dat;

    // inherited
    const bool& all_default_weights;
    const uint32_t verbosity;
};

struct PepinInt {
    PepinInt(const double _epsilon, const double _delta, const uint32_t seed,
              const uint32_t verbosity = 1);
    ~PepinInt();

    void set_force_eager(const int _force_eager) {
        force_eager = _force_eager;
    }
    void set_fast_center_calc(const int _fast_center_calc) {
        fast_center_calc = _fast_center_calc;
    }

    uint32_t new_vars(uint32_t n) {
        release_assert(nvars == 0);
        if ((n%4) != 0) {
            num_fake_vars = 4-(n%4);
            n += num_fake_vars;
        }
        print_verb(2, "fake num vars: " << num_fake_vars);

        nvars = n;
        weights.resize(nvars);
        bucket.set_nvars(nvars);
        seen.resize(nvars, false);

        return nvars;
    }

    bool add_clause(const vector<Lit>& cl);
    void magic(const vector<Lit>& cl, mpz_t ni);
    void get_cl_precision(const vector<Lit>& cl, mpz_t cl_prec_out);
    void poisson(mpz_t n_local, mpq_t sampl_prob, mpz_t samples_needed_out);
    void add_samples(const vector<Lit>& cl, const uint64_t dnf_cl_num, const uint64_t num_samples);

    void check_ready() const;
    const mpf_t* get_low_prec_appx_num_points() const;
    const mpf_t* get_low_prec_appx_weighted_sol() const;
    const mpq_t* get_appx_weighted_sol() const;

    uint32_t nVars() const {
        return nvars;
    }

    inline void set_var_weight(
            const uint32_t var,
            const uint32_t dividend,
            const uint32_t divisor)
    {
        if (!(dividend == 1 && divisor == 2)) all_default_weights = false;
        weights[var].dividend = dividend;
        weights[var].divisor = divisor;
        weights[var].bits_of_entropy =  weights[var].calc_bits_of_entropy();
        print_verb(2, "var " << (var+1)
        << " dividend:" << dividend
        << " divisor:" << divisor);
    }

    void set_n_cls(uint32_t n_cls);
    const char* get_version_info() const;
    const char* get_compilation_env() const;

    // For calculating birthday paradox
    vector<Weight> weights;
    bool all_default_weights = true;
    vector<uint8_t> seen;

    uint32_t num_fake_vars = 0; // needed to be able to have nvars %4 == 0
    uint32_t nvars = 0;
    uint32_t n_cls_declared = 0;
    double epsilon;
    double delta;
    uint64_t thresh = 0; //intentionally wrong
    mpq_t sampl_prob;
    uint32_t sampl_prob_expbit = 0;
    uint32_t sampl_prob_expbit_before_approx = std::numeric_limits<uint32_t>::max();
    uint32_t sampl_prob_expbit_before_magic = std::numeric_limits<uint32_t>::max();
    mpz_t prod_precision;
    Bucket bucket;
    uint64_t num_cl_added = 0;
    uint32_t verbosity;
    bool force_eager = false;
    double last_10k_time;
    int fast_center_calc = true;
    uint64_t added_samples_during_processing = 0;
    std::mt19937_64 mtrand;

    uint64_t lazy_samples_called = 0;
    uint64_t samples_called = 0;

    //Values that are re-calculated ONLY when needed
    mpz_t last_n_magic;
    mpq_t center_magic;
    mpz_t last_n_appx;
    mpq_t center_for_a;

    //constant values
    mpq_t constant_one;
    mpz_t constant_one_z;
    mpq_t constant_half;
    mpq_t constant_hundred;
    mpq_t constant_very_small_center;

    //temporaries
    mpq_t center_100;
    mpz_t n;
    mpz_t ni_plus_bucketsz;
    mpz_t ni;
    vector<Lit> cl_tmp;


    // Return
    bool ret_set = false;
    mpq_t weigh_num_sols;
    mpf_t low_prec_weigh_num_sols;
    mpf_t low_prec_num_points;
};

}

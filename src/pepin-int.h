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
#include <utility>
#include <algorithm>
#include <gmpxx.h>
#include <string.h>
#include "pepin.h"

using std::vector;
using std::pair;
using std::make_pair;
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

// Sample-indexed storage. Each "sample" holds nvars logical positions,
// each in {0, 1, 3=unset}. DenseElems packs all positions as 2-bit cells;
// SparseElems stores only positions that have been pinned to a concrete value.
class DenseElems {
public:
    void set_nvars(uint32_t _nvars) { nvars = _nvars; }
    uint64_t num_samples() const { return n_samples; }

    value get(uint64_t sample_idx, uint32_t var) const {
        const uint64_t at = sample_idx*(uint64_t)nvars + var;
        const uint64_t at_pos = at >> 2;
        const uint32_t sub_at = (at & 3) << 1;
        return (data[at_pos] >> sub_at) & 3;
    }

    void set(uint64_t sample_idx, uint32_t var, value val) {
        assert(val < 3);
        const uint64_t at = sample_idx*(uint64_t)nvars + var;
        const uint64_t at_pos = at >> 2;
        const uint32_t sub_at = (at & 3) << 1;
        data[at_pos] &= ~((uint8_t)3 << sub_at);
        data[at_pos] |= (uint8_t)(val << sub_at);
    }

    // Append a new sample (all positions unset). Amortized O(nvars/4)
    // thanks to vector's exponential growth.
    uint64_t add_sample() {
        assert(nvars % 4 == 0);
        const uint64_t old_size = data.size();
        data.resize(old_size + nvars/4, 0xff);
        return n_samples++;
    }

    // Reset a previously allocated sample slot to all-unset (slot reuse).
    void reset_sample(uint64_t sample_idx) {
        assert(nvars % 4 == 0);
        const uint64_t at_pos = sample_idx*(uint64_t)nvars/4;
        memset(data.data()+at_pos, 0xff, nvars/4);
    }

private:
    uint32_t nvars = 0;
    uint64_t n_samples = 0;
    std::vector<uint8_t> data;
};

// Hybrid per-sample storage. Each sample starts in sparse mode (sorted
// vector of (var, value) pairs, binary-searched). Once its concrete-entry
// count crosses nvars / DENSE_SWITCH_RATIO, it auto-promotes to a packed
// 2-bit-per-var dense bitset, at which point get/set become O(1).
//
// Short-lived samples (killed quickly by remove_sol) never promote and
// stay cheap. Long-lived samples that would otherwise accumulate O(N^2)
// insert cost in a sorted vector pay a one-shot promotion and then run
// at dense speed.
class SparseElems {
    static bool cmp_var(const std::pair<uint32_t, uint8_t>& p, uint32_t v) {
        return p.first < v;
    }

    struct Sample {
        std::vector<std::pair<uint32_t, uint8_t>> sparse;
        std::vector<uint8_t> dense; // empty until promoted; size nvars/4
        bool is_dense = false;
    };

    // Promote when sparse entries exceed this. Tuned for runtime, not
    // memory: sparse insert is O(N) per shift, so total per-sample insert
    // cost is O(N^2); dense set is O(1). Crossover where per-insert cost
    // of sparse exceeds dense is roughly N>=8, so we cap N around there.
    // For very long-lived samples this means we briefly stay sparse, then
    // run at dense speed for the bulk of remove_sol passes. Short-lived
    // samples (high k -> low survival, or k=3 with rapid kills) never
    // promote and keep the sparse-side wins of cheap add/reset.
    static constexpr uint32_t DENSE_SWITCH_THRESH = 32;

    void promote(Sample& s) const {
        s.dense.assign(nvars/4, 0xff);
        for (const auto& p : s.sparse) {
            const uint64_t at_pos = p.first >> 2;
            const uint32_t sub_at = (p.first & 3) << 1;
            s.dense[at_pos] &= ~((uint8_t)3 << sub_at);
            s.dense[at_pos] |= (uint8_t)(p.second << sub_at);
        }
        std::vector<std::pair<uint32_t, uint8_t>>().swap(s.sparse);
        s.is_dense = true;
    }

public:
    void set_nvars(uint32_t _nvars) {
        assert(_nvars % 4 == 0);
        nvars = _nvars;
    }
    uint64_t num_samples() const { return samples.size(); }

    value get(uint64_t sample_idx, uint32_t var) const {
        const Sample& s = samples[sample_idx];
        if (s.is_dense) {
            const uint64_t at_pos = var >> 2;
            const uint32_t sub_at = (var & 3) << 1;
            return (s.dense[at_pos] >> sub_at) & 3;
        }
        auto it = std::lower_bound(s.sparse.begin(), s.sparse.end(), var, cmp_var);
        if (it != s.sparse.end() && it->first == var) return it->second;
        return 3;
    }

    void set(uint64_t sample_idx, uint32_t var, value val) {
        assert(val < 3);
        Sample& s = samples[sample_idx];
        if (s.is_dense) {
            const uint64_t at_pos = var >> 2;
            const uint32_t sub_at = (var & 3) << 1;
            s.dense[at_pos] &= ~((uint8_t)3 << sub_at);
            s.dense[at_pos] |= (uint8_t)(val << sub_at);
            return;
        }
        auto it = std::lower_bound(s.sparse.begin(), s.sparse.end(), var, cmp_var);
        if (it != s.sparse.end() && it->first == var) {
            it->second = (uint8_t)val;
            return;
        }
        s.sparse.emplace(it, var, (uint8_t)val);
        if (s.sparse.size() > DENSE_SWITCH_THRESH) promote(s);
    }

    uint64_t add_sample() {
        samples.emplace_back();
        return samples.size() - 1;
    }

    void reset_sample(uint64_t sample_idx) {
        Sample& s = samples[sample_idx];
        // Reset to sparse mode; release the dense buffer if any so the
        // freed slot doesn't carry stale promotion state into reuse.
        s.sparse.clear();
        std::vector<uint8_t>().swap(s.dense);
        s.is_dense = false;
    }

private:
    uint32_t nvars = 0;
    std::vector<Sample> samples;
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

template<typename T>
class Bucket {
public:
    Bucket (const bool& _all_default_weights, const uint32_t& _verbosity) :
        all_default_weights(_all_default_weights),
        verbosity(_verbosity)
    {}

    // Returns a sample_idx (not a linear position).
    uint64_t add_lazy_common(const uint64_t dnf_cl_num) {
        uint64_t at;
        if (empties.empty()) {
            at = elems.add_sample();
            elems_dat.push_back(ElemDat(dnf_cl_num));
        } else {
            at = empties.back();
            empties.pop_back();
            elems.reset_sample(at);
            elems_dat[at] = ElemDat(dnf_cl_num);
        }
        return at;
    }

    void add_lazy(const vector<Lit>& cl, const uint64_t dnf_cl_num);
    uint64_t get_size() const { return size; }
    size_t nVars() const { return nvars; }

    void remove_sol(const vector<Lit>& cl,
                    const vector<Weight>& var_weights,
                    std::mt19937_64& mtrand);

    void remove(const size_t at) {
        assert(elems_dat[at].empty == false);
        assert(elems.num_samples() > at);

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
        elems.set_nvars(_nvars);
    }

    void print_elems_stats(const uint64_t tot_num_dnf_cls) const;

private:
    uint64_t size = 0;
    uint32_t nvars = 0;
    T elems;
    vector<uint32_t> empties;
    vector<ElemDat> elems_dat;

    // inherited
    const bool& all_default_weights;
    const uint32_t verbosity;
};

// Base class for runtime polymorphism
struct PepinIntBase {
    virtual ~PepinIntBase() {}
    virtual void set_force_eager(const int _force_eager) = 0;
    virtual void set_fast_center_calc(const int _fast_center_calc) = 0;
    virtual uint32_t new_vars(uint32_t n) = 0;
    virtual bool add_clause(const vector<Lit>& cl) = 0;
    virtual uint32_t nVars() const = 0;
    virtual void set_var_weight(const uint32_t var, const uint32_t dividend, const uint32_t divisor) = 0;
    virtual void set_n_cls(uint32_t n_cls) = 0;
    virtual const mpf_t* get_low_prec_appx_num_points() const = 0;
    virtual const mpf_t* get_low_prec_appx_weighted_sol() const = 0;
    virtual const mpq_t* get_appx_weighted_sol() const = 0;
};

template<typename StorageType>
class PepinInt : public PepinIntBase {
public:
    PepinInt(const double _epsilon, const double _delta, const uint32_t seed,
              const uint32_t verbosity = 1);
    ~PepinInt();

    void set_force_eager(const int _force_eager) override {
        force_eager = _force_eager;
    }
    void set_fast_center_calc(const int _fast_center_calc) override {
        fast_center_calc = _fast_center_calc;
    }

    uint32_t new_vars(uint32_t n) override {
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

    bool add_clause(const vector<Lit>& cl) override;
    void magic(const vector<Lit>& cl, mpz_t ni);
    void get_cl_precision(const vector<Lit>& cl, mpz_t cl_prec_out);
    void poisson(mpz_t n_local, mpq_t sampl_prob, mpz_t samples_needed_out);
    void add_samples(const vector<Lit>& cl, const uint64_t dnf_cl_num, const uint64_t num_samples);

    void check_ready() const;
    const mpf_t* get_low_prec_appx_num_points() const override;
    const mpf_t* get_low_prec_appx_weighted_sol() const override;
    const mpq_t* get_appx_weighted_sol() const override;

    uint32_t nVars() const override { return nvars; }
    inline void set_var_weight(
            const uint32_t var,
            const uint32_t dividend,
            const uint32_t divisor) override
    {
        if (!(dividend == 1 && divisor == 2)) all_default_weights = false;
        weights[var].dividend = dividend;
        weights[var].divisor = divisor;
        weights[var].bits_of_entropy =  weights[var].calc_bits_of_entropy();
        print_verb(2, "var " << (var+1)
        << " dividend:" << dividend
        << " divisor:" << divisor);
    }

    void set_n_cls(uint32_t n_cls) override;
    const char* get_version_info() const;
    const char* get_compilation_env() const;

private:

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
    Bucket<StorageType> bucket;
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

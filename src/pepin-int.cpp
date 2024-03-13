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

#include "GitSHA1.h"
#include "pepin-int.h"
#include "time_mem.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <cmath>
#include <string.h>

using std::cout;
using std::endl;

using namespace PepinIntNS;

void Bucket::remove_half(std::mt19937_64& mtrand)
{
    const auto orig_size = get_size();
    if (get_size() > 0) {
        print_verb(1, "Removing half ... bucket size now: " << get_size());
    }
    for(size_t i = 0, at = 0; i < elems.size(); i+=nvars, at++) {
        if (elems_dat[at].empty) continue;
        int ok = mtrand()&1;
        if (ok) continue;
        remove(at);
    }
    if (orig_size > 0) {
        print_verb(2, " after: " << get_size() << " ratio: " << std::setprecision(2) << ((double)get_size()/(double)orig_size));
    }
}

void Bucket::remove_sol(const vector<Lit>& cl,
                        const vector<Weight>& var_weights,
                        std::mt19937_64& mtrand)
{
    const size_t size_before = get_size();
    uint64_t rand_pool;
    uint32_t rand_pool_bits = 0;
    for(size_t i = 0, at = 0; i < elems.size(); i+=nvars, at++) {
        if (elems_dat[at].empty) continue;
        elems_dat[at].checked++;

        bool sol = true;
        bool updated = false;
        for(const auto& lit: cl) {
            const uint32_t var = lit.var();

            //lazy, let's make it real
            if (elems[i+var] == 3) {
                updated = true;
                WeightPrec w;

                if (all_default_weights || var_weights[var].divisor == 2) {
                    if (rand_pool_bits == 0) {rand_pool = mtrand(); rand_pool_bits = 64;}
                    w = rand_pool & 1;
                    rand_pool_bits--;
                    rand_pool >>= 1;
                } else {
                    std::uniform_int_distribution<int> uid(0,var_weights[var].divisor-1);
                    w = uid(mtrand);
                }
                if ((all_default_weights && (w == 0))
                    || (!all_default_weights && w < var_weights[var].dividend)) {
                    elems.set(i+var, 1);
                    if (!lit.sign() != 1) {sol=false;break;}
                } else {
                    elems.set(i+var, 0);
                    if (!lit.sign() != 0) {sol=false;break;}
                }
            } else {
                //Filter on eager
                if (elems[i+var] != !lit.sign()) {sol=false;break;}
            }
        }
        if (sol) remove(at);
        else elems_dat[at].updated+=updated;
    }
    if (size_before > get_size()) {
        print_verb(2, "Filtered. Removed: " << size_before - get_size());
    }
}

void Bucket::print_contents() const
{
    cout << "-- Bucket contents -- SZ: " << get_size() << endl;
    for(uint32_t i = 0, at = 0; i < elems.size(); i+=nvars, at++) {
        if (elems_dat[at].empty) continue;
        cout << "[" << std::setw(5) << at << "]: ";
        for(uint32_t i2 = 0; i2 < nvars; i2++) {
            cout << (int)elems[i+i2];
        }
        cout << endl;
    }
    cout << "-- Bucket contents end --: " << endl;
}

void PepinInt::set_n_cls(uint32_t n_cls)
{
    thresh = (int)(4.0*(std::log2(n_cls+1)/
        (epsilon*epsilon))*std::log2(1.0/delta));
    print_verb(2, "Threshold computed: " << thresh
    << " -- cl num: " << n_cls
    << " epsilon: " << std::setprecision(5) << epsilon
    << " delta: " << std::setprecision(5) <<delta);

    n_cls_declared = n_cls;
}

PepinInt::PepinInt(const double _epsilon, const double _delta,
                     const uint32_t seed, const uint32_t _verbosity) :
    bucket(all_default_weights, _verbosity)
{
    mtrand.seed(seed);
    epsilon = _epsilon;
    delta = _delta;
    verbosity = _verbosity;

    // set sample probability to 1
    mpq_init(sampl_prob);

    //set prod_precision to 1
    mpz_init(prod_precision);
    mpz_set_ui(prod_precision, 1);

    //init ni
    mpz_init(ni);

    mpq_init(constant_one);
    mpq_set_ui(constant_one, 1, 1);

    mpq_init(center_for_a);
    mpq_init(center_magic);
    mpq_init(constant_half);
    mpq_set_ui(constant_half, 1, 2);
    mpz_init(n);
    mpz_init(last_n_appx);
    mpz_init(last_n_magic);

    mpq_init(constant_hundred);
    mpq_set_ui(constant_hundred, 100, 1);

    mpq_init(constant_very_small_center);
    mpq_set_ui(constant_very_small_center, 1, 1UL<<31);

    mpq_init(center_100);
    mpz_init(ni_plus_bucketsz);

    mpz_init(constant_one_z);
    mpz_set_ui(constant_one_z, 1);
}

PepinInt::~PepinInt()
{
    mpq_clear(constant_one);
    mpz_clear(ni);
    mpz_clear(prod_precision);
    mpq_clear(sampl_prob);
    mpq_clear(center_for_a);
    mpq_clear(center_magic);
    mpq_clear(constant_very_small_center);
    mpq_clear(constant_half);
    mpz_clear(n);
    mpz_clear(last_n_appx);
    mpz_clear(last_n_magic);
    mpq_clear(constant_hundred);
    mpq_clear(center_100);
    mpz_clear(ni_plus_bucketsz);
    mpz_clear(constant_one_z);

    // Return values freed
    mpq_clear(weigh_num_sols);
    mpf_clear(low_prec_weigh_num_sols);
    mpf_clear(low_prec_num_points);
}

const char* PepinInt::get_version_info() const {
    return PepinIntNS::get_version_sha1();
}

const char* PepinInt::get_compilation_env() const {
    return PepinIntNS::get_compilation_env();
}

void PepinInt::get_cl_precision(const vector<Lit>& cl, mpz_t cl_prec_out)
{
    bool fast_ok = true;
    for(const Lit& l: cl) {
        if (weights[l.var()].divisor == 2 &&
            weights[l.var()].dividend == 1)
        {
            //
        } else {
            fast_ok = false;
            break;
        }
    }
    if (fast_ok) {
        mpz_mul_2exp(cl_prec_out, constant_one_z, nvars-cl.size());
    } else {
        assert(cl.size() > 0);
        mpz_set(cl_prec_out, prod_precision);

        mpz_t x;
        mpz_init(x);
        mpz_set_ui(x, weights[cl[0].var()].divisor);
        for(uint32_t i = 1; i < cl.size(); i++) {
            uint32_t var = cl[i].var();
            mpz_mul_ui(x, x, weights[var].divisor);
        }
        mpz_divexact(cl_prec_out, cl_prec_out, x);
        mpz_clear(x);

        for(const Lit l: cl) {
            uint32_t var = l.var();
            if (!l.sign()) {
                mpz_mul_ui(cl_prec_out, cl_prec_out, weights[var].dividend);
            } else {
                mpz_mul_ui(cl_prec_out, cl_prec_out, weights[var].divisor-weights[var].dividend);
            }
        }
    }
}

void PepinInt::magic(const vector<Lit>& cl, mpz_t samples_needed)
{
    get_cl_precision(cl, n);

    print_verb(3, "Clause precision: " << n);

    //calculate center_magic = n*sampl_prob
    if (sampl_prob_expbit_before_magic == sampl_prob_expbit &&
        mpz_cmp(last_n_magic, n) == 0 &&
        fast_center_calc)
    {
        //nothing to do, reuse old center_magic and sampl_prob values
    } else {
        mpq_div_2exp(sampl_prob, constant_one, sampl_prob_expbit); ///recreate sampl_prob
        mpq_set_z(center_magic, n);
        mpq_mul(center_magic, center_magic, sampl_prob);
        sampl_prob_expbit_before_magic = sampl_prob_expbit;
        mpz_set(last_n_magic, n);
    }

    // Half things until it's safe to do so
    bool quick_halfing_ran = false;
    while (mpq_cmp_ui(center_magic, thresh, 1) > 0) {
        //mpq_cmp_ui(u,vn,vd) -- Compare U with Vn/Vd
        quick_halfing_ran = true;
        print_verb(3, "center_magic is larger than thresh: "
                << thresh << " halfing, incrementing sampl_prob_expbit.");
        mpq_mul(center_magic, center_magic, constant_half);
        sampl_prob_expbit++;
        print_verb(2, "Sampl. prob is now: 2**-" << sampl_prob_expbit
            << " nvars-expprob = " << nvars-sampl_prob_expbit
            << " bucket size: " << bucket.get_size());
        bucket.remove_half(mtrand);
    }
    if (quick_halfing_ran) {
        print_verb(2, "Sampl. prob is now: 2**-" << sampl_prob_expbit
        << " nvars-expprob = " << nvars-sampl_prob_expbit
        << " bucket size: " << bucket.get_size());
    }
    print_verb(2, "Approx is: " << center_magic);
    mpz_set_ui(ni_plus_bucketsz, 0);

    ///recreate sampl_prob
    mpq_div_2exp(sampl_prob, constant_one, sampl_prob_expbit);
    approx_binomial(n, sampl_prob, samples_needed);

    while(true)  {
        mpz_add_ui(ni_plus_bucketsz, samples_needed, bucket.get_size());
        print_verb(2, "bucket size: " << bucket.get_size());
        print_verb(2, "ni_plus_bucketsz:" << ni_plus_bucketsz);
        if (mpz_cmp_ui(ni_plus_bucketsz, thresh) < 0) {
            print_verb(3, "Less than threshold " << thresh << " -> breaking");
            break;
        }

        approx_binomial(samples_needed, constant_half, samples_needed);
        sampl_prob_expbit++;
        print_verb(2, "Sampl. prob is now: 2**-" << sampl_prob_expbit
        << " nvars-expprob = " << nvars-sampl_prob_expbit
        << " bucket size: " << bucket.get_size());
        bucket.remove_half(mtrand);
    }
}

void PepinInt::approx_binomial(
    mpz_t n_local, mpq_t sampl_prob, mpz_t samples_needed_out)
{
    
    mpq_set_z(center_for_a, n_local);
    mpq_mul(center_for_a, center_for_a, sampl_prob);
    const double center_low_prec = mpq_get_d(center_for_a);
    assert(center_low_prec != nan(""));
    assert(!isinf(center_low_prec));

    std::poisson_distribution<> poiss(center_low_prec);
    uint64_t ni = poiss(mtrand);
    mpz_set_ui(samples_needed_out, ni);

    print_verb(2, "Used poission distribution");
    print_verb(2, "Mu (i.e. center): " << center_low_prec);
    return;
    
}

struct SampleSorter {
    SampleSorter(const uint32_t _nvars, const vector<value>& _sol,
                 const vector<WeightPrec>& _w):
                 w(_w), sol(_sol), nvars(_nvars)
    {
    }

    bool operator()(const Sample& s1, const Sample& s2) const {
        for(uint32_t i = 0; i < nvars; i ++) {
            if (sol[s1.sol_at+i] != sol[s2.sol_at+i]) {
                return sol[s1.sol_at+i] < sol[s2.sol_at+i];
            }
        }
        for(uint32_t i = 0; i < nvars; i ++) {
            if (w[s1.ws_at+i] != w[s2.ws_at+i]) {
                return w[s1.ws_at+i] < w[s2.ws_at+i];
            }
        }

        assert(s1.gen_point != s2.gen_point);
        return s1.gen_point < s2.gen_point;
    }

    const vector<WeightPrec>& w;
    const vector<value>& sol;
    const uint32_t nvars;
};

struct SampleSorterGenpoint {
    bool operator()(const Sample& s1, const Sample& s2) const {
        //return TRUE when s1.gen_point is earlier
        return s1.gen_point < s2.gen_point;
    }
};

void Bucket::add_lazy(const vector<Lit>& cl, const uint64_t dnf_cl_num)
{
    const uint64_t at = add_lazy_common(dnf_cl_num);
    for(const Lit l: cl) elems.set(at+l.var(), !l.sign());
    size++;
}

void PepinInt::add_uniq_samples(
        const vector<Lit>& cl, const uint64_t dnf_cl_num, const uint64_t num)
{
    samples_called++;
    double bits_of_entropy = 0;
    if (all_default_weights) {
        bits_of_entropy = ((double)nvars-(double)cl.size());
    } else {
        for(const auto& l: cl) seen[l.var()] = true;
        for(uint32_t i = 0; i < nVars(); i++) {
            if (seen[i]) continue;
            bits_of_entropy += weights[i].bits_of_entropy;
        }
        for(const auto& l: cl) seen[l.var()] = false;
    }
    //Even accounting for birthday paradox, there are still 30 bits left
    bool lazy = bits_of_entropy/2.0-std::log2(num+1) > 30;

    if (force_eager) lazy = false;
    print_verb(3, "Lazy: " << lazy);
    lazy_samples_called += lazy;

    if (lazy) {
        added_samples_during_processing += num;
        for(uint64_t i = 0; i < num; i ++) {
            bucket.add_lazy(cl, dnf_cl_num);
        }
        return;
    }

    //We'll add the samples into these, samples[] will be the pointer
    vector<value> sol;
    vector<WeightPrec> ws;
    vector<Sample> samples;

    double myTime = cpuTime();
    uint64_t todo = num + 1;
    uint64_t rand_pool;
    uint32_t rand_pool_bits = 0;
    uint64_t at = 0;

    uint32_t retries = 0;
    while (samples.size() < num) {
        if (retries > 100) {
            cout << "ERROR. This version of the algorithm can only deal with unique samples, and the precision requested would require more samples than there is volume. So we can't do that. Please lower your epsilon." << endl;
            exit(-1);
        }
        sol.resize(sol.size()+nVars()*todo);
        ws.resize(ws.size()+nVars()*todo);
        for(uint64_t i = 0; i < todo; i ++) {
            Sample s;
            size_t start_at = at;

            //Generate completely random values, then fix later to
            //set the clause's values (that are very few)
            //  --> faster this way, less branching here
            for(uint32_t var = 0; var < nvars; var ++) {
                //Set up the sample data
                if (var == 0) {
                    s.sol_at = at/nvars;
                    s.ws_at = at/nvars;
                    s.gen_point = samples.size();
                }

                WeightPrec w;
                if (weights[var].divisor == 2 &&
                    weights[var].dividend == 1)
                {
                    //use fast bit system if possible
                    if (rand_pool_bits == 0) {
                        rand_pool = mtrand();
                        rand_pool_bits = 64;
                    }
                    w = rand_pool & 1;
                    rand_pool_bits--;
                    rand_pool >>= 1;
                } else {
                    std::uniform_int_distribution<int> uid(0,weights[var].divisor-1);
                    w = uid(mtrand);

                }
                sol[at] = w < weights[var].dividend;
                ws[at] = w;
                at++;
            }

            //Fix the clause values now
            for(const Lit lit: cl) {
                uint32_t var = lit.var();
                size_t my_at = start_at + var;
                sol[my_at] = !lit.sign();
                if (sol[my_at]) {
                    std::uniform_int_distribution<int> uid(0,weights[var].dividend-1);
                    ws[my_at] = uid(mtrand);
                } else {
                    std::uniform_int_distribution<int> uid(0,weights[var].divisor-weights[var].dividend-1);
                    ws[my_at] = weights[var].dividend+uid(mtrand);
                }
            }

            samples.push_back(s);
        }
        print_verb(2, "Generated " << todo << " extra non-unique samples");
        assert(sol.size() == ws.size());
        assert(sol.size() % nVars() == 0);

        //remove duplicates
        if (!samples.empty()) {
            std::sort(samples.begin(), samples.end(), SampleSorter(nvars, sol, ws));
            uint64_t j = 0;
            for(uint64_t i = 1; i < samples.size(); i++) {
                if (!samples[i].equals_sol_w(samples[j], nvars, sol, ws)) {
                    j++;
                    samples[j] = samples[i];
                }
            }
            samples.resize(j+1);
        }
        if (samples.size() < num) todo = (num - samples.size())*2;
        print_verb(2, "Now we have " << samples.size() << " unique samples, T:"
            << (cpuTime() - myTime));
        retries++;
    }

    //add unique samples
    assert(samples.size() >= num);
    if (!lazy) {
        std::sort(samples.begin(), samples.end(), SampleSorterGenpoint());
    }
    for(uint64_t i = 0; i < num; i ++) {
        bucket.add(&sol[samples[i].sol_at], dnf_cl_num);
    }

    print_verb(2, "Added " << num << " unique samples, new bucket size: " << bucket.get_size());
}

void Bucket::print_elems_stats(const uint64_t tot_num_dnf_cls) const
{
    uint64_t num_elems = 0;
    uint64_t sum_dnf_cl_num = 0;
    uint64_t sum_updated = 0;
    uint64_t sum_checked = 0;
    uint64_t min_dnf_cl_num = std::numeric_limits<uint64_t>::max();
    uint64_t max_dnf_cl_num = 0;
    for(const auto& dat: elems_dat) {
        if (dat.empty) continue;
        num_elems++;
        sum_dnf_cl_num += dat.dnf_cl_num;
        sum_updated += dat.updated;
        sum_checked += dat.checked;
        min_dnf_cl_num = std::min(min_dnf_cl_num, dat.dnf_cl_num);
        max_dnf_cl_num = std::max(max_dnf_cl_num, dat.dnf_cl_num);
    }
    cout << "Num elems: " << num_elems
    << " avg dnf_cl_num: " << std::setprecision(2) << std::fixed << (double)sum_dnf_cl_num/(double)num_elems
    << " avg num updated: " << std::setprecision(2) << std::fixed << (double)sum_updated/(double)num_elems
    << " avg num checked: " << std::setprecision(2) << std::fixed << (double)sum_checked/(double)num_elems
    << " max age: " << (tot_num_dnf_cls-min_dnf_cl_num)
    << " min age: " << (tot_num_dnf_cls-max_dnf_cl_num)
    << endl;
}

bool PepinInt::add_clause(const vector<Lit>& cl) {
    assert(thresh != 0 && "The number of clauses was not set beforehand!");
    assert(num_cl_added < n_cls_declared);
    if (num_cl_added == 0) {
        print_verb(2, "First clause added, calculating prod precision");
        for(const auto w: weights) {
            mpz_mul_ui(prod_precision, prod_precision, w.divisor);
        }
        last_10k_time = cpuTime();
    }
    if (num_cl_added >= n_cls_declared) {
        cout << "ERROR: you are trying to add more clauses than were declared." << endl;
        assert(false);
        exit(-1);
    }

    print_verb(2, "Adding clause: " << cl);
    print_verb(2, "CL num: " << num_cl_added);
    assert(bucket.nVars() == nVars());
    for(const Lit& l: cl) assert(l.var() < nvars);

    // No duplicates
    cl_tmp = cl;
    sort(cl_tmp.begin(), cl_tmp.end());
    cl_tmp.erase(unique( cl_tmp.begin(), cl_tmp.end() ), cl_tmp.end());

    print_verb(2, "Filtering bucket from solutions. Orig sz:" << bucket.get_size());
    bucket.remove_sol(cl_tmp, weights, mtrand);
    print_verb(2, "Bucket size now: " << bucket.get_size());
    if (verbosity >= 3) bucket.print_contents();

    // Computing number of samples needed
    magic(cl_tmp, ni);

    print_verb(2, "Generating " << ni << " samples");
    if (!mpz_fits_ulong_p(ni)) {
        cout << "ERROR: ni does not fit a long unsigned int! Value is:" << endl;
        mpz_out_str(stdout, 10, ni);
        cout << endl;
        exit(-1);
    }
    uint64_t num_samples = mpz_get_ui(ni);
    if (num_samples > 0) add_uniq_samples(cl_tmp, num_cl_added, num_samples);

    num_cl_added++;
    if (verbosity >= 2) {
        cout << "-- after add_clause --" << endl;
        mpq_div_2exp(sampl_prob, constant_one, sampl_prob_expbit);
        //cout << "sampl_prob:" << sampl_prob << endl;
        cout << "bucket size:" << bucket.get_size() << endl;
        if (verbosity >= 3) bucket.print_contents();
        cout << "Threshold: " << thresh << endl;
    }
    if ((num_cl_added & 0x3fff) == 0x3fff) {
        print_verb(1, "--==>> Num CL processed: " << std::setw(10) << num_cl_added
        << " took T:" << std::setw(7) << std::setprecision(4) << cpuTime()-last_10k_time
        << " Bucket sz: " << std::setw(5) << bucket.get_size()
        << " sampl added since last print: " << added_samples_during_processing);
        last_10k_time = cpuTime();
        added_samples_during_processing = 0;
    }

    if (num_cl_added == n_cls_declared) {
        print_verb(1, "*** Finished.***");
        print_verb(1, "Num CL processed: " << num_cl_added);
        print_verb(1, "Sample function called: " << samples_called
        << " of which lazy: " << std::setprecision(6)
        << ((double)lazy_samples_called/(double)samples_called)*100.0
        << " %");

        //Statistics about the samples
        if (verbosity >= 1) bucket.print_elems_stats(num_cl_added);

        //Exp description
        double buck_size_log2 = log2(bucket.get_size());
        buck_size_log2 += sampl_prob_expbit;
        print_verb(1, "Bucket expbit: " << sampl_prob_expbit << " bucket log2: " << log2(bucket.get_size()));
        print_verb(1, "bucket_size/sampl_prob : 2**" << buck_size_log2);

        // Calculate the number of points in the N-dimensional space
        // NOTE: Bignum cannot do 2**(float), so we do:
        //    z*2**floor(float)
        //    where we calculate z = 2**(float-floor(float)) in C++, i.e. low precision
        mpf_init2(low_prec_num_points, 1000);
        mpf_set_d(low_prec_num_points, std::exp2(buck_size_log2-std::floor(buck_size_log2)));
        mpf_mul_2exp(low_prec_num_points, low_prec_num_points, std::floor(buck_size_log2)-num_fake_vars);

        //bucket_size/sampl_prob
        mpq_init(weigh_num_sols);
        mpq_set_ui(weigh_num_sols, bucket.get_size(), 1);
        mpq_div_2exp(sampl_prob, constant_one, sampl_prob_expbit);
        mpq_div(weigh_num_sols, weigh_num_sols, sampl_prob);

        // (bucket_size/sampl_prob)/prod_precision
        mpq_t tmp2;
        mpq_init(tmp2);
        mpq_set_z(tmp2, prod_precision);
        mpq_div(weigh_num_sols, weigh_num_sols, tmp2);
        mpq_clear(tmp2);

        mpf_init2(low_prec_weigh_num_sols, 100);
        mpf_set_q(low_prec_weigh_num_sols, weigh_num_sols);

        ret_set = true;
        return true;
    }
    return false;
}

void PepinInt::check_ready() const
{
    if (!ret_set) {
        cout << "ERROR, return value not yet available" << endl;
        assert(false);
        exit(-1);
    }
}

const mpf_t* PepinInt::get_low_prec_appx_num_points() const {
    check_ready();
    return &low_prec_num_points;
}
const mpf_t* PepinInt::get_low_prec_appx_weighted_sol() const {
    check_ready();
    return &low_prec_weigh_num_sols;
}
const mpq_t* PepinInt::get_appx_weighted_sol() const {
    check_ready();
    return &weigh_num_sols;
}

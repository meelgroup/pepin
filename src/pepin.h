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

#include <cstdint>
#include <cassert>
#include <vector>
#include <ostream>
#include <gmpxx.h>

namespace PepinNS {
struct PepinPrivate;

constexpr uint32_t var_Undef(0xffffffffU >> 4);

class Lit
{
    uint32_t x;
    constexpr explicit Lit(uint32_t i) : x(i) { }
public:
    constexpr Lit() : x(var_Undef<<1) {}   // (lit_Undef)
    constexpr explicit Lit(uint32_t var, bool is_inverted) :
        x(var + var + is_inverted)
    {}

    constexpr Lit  operator~() const {
        return Lit(x ^ 1);
    }
    constexpr Lit  operator^(const bool b) const {
        return Lit(x ^ (uint32_t)b);
    }
    Lit& operator^=(const bool b) {
        x ^= (uint32_t)b;
        return *this;
    }
    constexpr bool sign() const {
        return x & 1;
    }
    constexpr uint32_t  var() const {
        return x >> 1;
    }
    constexpr Lit  unsign() const {
        return Lit(x & ~1U);
    }
    constexpr bool operator==(const Lit& p) const {
        return x == p.x;
    }
    constexpr bool operator!= (const Lit& p) const {
        return x != p.x;
    }
    /**
    @brief ONLY to be used for ordering such as: a, b, ~b, etc.
    */
    constexpr bool operator <  (const Lit& p) const {
        return x < p.x;     // '<' guarantees that p, ~p are adjacent in the ordering.
    }
    constexpr bool operator >  (const Lit& p) const {
        return x > p.x;
    }
    constexpr bool operator >=  (const Lit& p) const {
        return x >= p.x;
    }
};

constexpr Lit itol(int32_t data)
{
    assert(data != 0);
    return Lit(abs(data)-1, data < 0);
}

static const Lit lit_Undef(var_Undef, false);  // Useful special constants.

inline std::ostream& operator<<(std::ostream& os, const Lit lit)
{
    if (lit == lit_Undef) {
        os << "lit_Undef";
    } else {
        os << (lit.sign() ? "-" : "") << (lit.var() + 1);
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& co, const std::vector<Lit>& lits)
{
    for (uint32_t i = 0; i < lits.size(); i++) {
        co << lits[i];

        if (i != lits.size()-1)
            co << " ";
    }

    return co;
}

struct Pepin {
    Pepin(const double epsilon = 0.5, const double delta = 0.5, const uint32_t seed = 1,
              const uint32_t verbosity = 0);
    ~Pepin();

    void set_force_eager(const int force_eager);
    void set_fast_center_calc(const int fast_center_calc);
    uint32_t new_vars(const uint32_t n);
    bool add_clause(const std::vector<Lit>& cl);
    uint32_t nVars() const;
    void set_var_weight(
            const uint32_t var,
            const uint32_t dividend,
            const uint32_t divisor);

    void set_n_cls(uint32_t n_cls);
    static const char* get_version_sha1();
    static const char* get_compilation_env();
    const mpf_t* get_low_prec_appx_num_points() const;
    const mpf_t* get_low_prec_appx_weighted_sol() const;
    const mpq_t* get_appx_weighted_sol() const;

    PepinPrivate* pepin = NULL;
};

}

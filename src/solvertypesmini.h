/******************************************
Pepin

Copyright (C) 2021 Mate Soos, Kuldeep S. Meel

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
***********************************************/

#pragma once

#include <cstdint>
#include <iostream>
#include <cassert>
#include <vector>
#include <array>
#include <algorithm>

namespace PepinNS {

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
    constexpr static Lit toLit(uint32_t data)
    {
        return Lit(data);
    }
};

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

}

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
#include <solvertypesmini.h>
#include <vector>

namespace PepinNS {
struct PepinPrivate;

struct Pepin {
    Pepin(const double epsilon, const double delta, const uint32_t seed,
              const uint32_t verbosity = 1);
    ~Pepin();

    void set_force_eager(const int _force_eager);
    void set_fast(const int _fast);
    uint32_t new_vars(const uint32_t n);
    void add_clause(const std::vector<Lit>& cl, const uint64_t dnf_cl_num);
    void add_uniq_samples(const std::vector<Lit>& cl, const uint64_t dnf_cl_num, const uint64_t num_samples);
    uint32_t nVars() const;
    void set_var_weight(
            const uint32_t var,
            const uint32_t dividend,
            const uint32_t divisor);

    void set_n_cls(uint32_t n_cls);
    const char* get_version_info() const;
    const char* get_compilation_env() const;
    PepinPrivate* pepin;
};

}

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

#include "pepin.h"
#include "pepin-int.h"
#include "GitSHA1.h"

using namespace PepinIntNS;
using PepinNS::Lit;

struct PepinNS::PepinPrivate {
    PepinPrivate(PepinIntNS::PepinInt *_p) : p(_p) {}
    ~PepinPrivate() {delete p;}
    PepinIntNS::PepinInt* p = NULL;
};

PepinNS::Pepin::Pepin(const double epsilon, const double delta, const uint32_t seed,
          const uint32_t verbosity) {
    pepin = new PepinPrivate(new PepinIntNS::PepinInt(epsilon, delta, seed, verbosity));
}
PepinNS::Pepin::~Pepin() { delete pepin; }

void PepinNS::Pepin::set_force_eager(const int _force_eager) {
    pepin->p->set_force_eager(_force_eager);
}

void PepinNS::Pepin::set_fast(const int _fast) { pepin->p->set_fast(_fast); }
uint32_t PepinNS::Pepin::nVars() const { return pepin->p->nVars(); }
uint32_t PepinNS::Pepin::new_vars(const uint32_t n) { return pepin->p->new_vars(n); }

void PepinNS::Pepin::add_clause(const std::vector<Lit>& cl, const uint64_t dnf_cl_num)
{
    pepin->p->add_clause(cl, dnf_cl_num);
}

void PepinNS::Pepin::set_var_weight(
        const uint32_t var,
        const uint32_t dividend,
        const uint32_t divisor)
{
    pepin->p->set_var_weight(var, dividend, divisor);
}

void PepinNS::Pepin::set_n_cls(uint32_t n_cls) { pepin->p->set_n_cls(n_cls); }
const char* PepinNS::Pepin::get_version_info() { return PepinIntNS::get_version_sha1(); }
const char* PepinNS::Pepin::get_compilation_env() { return PepinIntNS::get_compilation_env(); }

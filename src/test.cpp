/*
 Pepin

 Copyright (c) 2021 All rights reserved. Authors:
    Mate Soos
    Kuldeep S. Meel

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
#include <iostream>
#include <iomanip>
#include <vector>
using namespace PepinNS;

int main() {
  Pepin pepin;
  pepin.new_vars(20);
  pepin.set_n_cls(2);

  typedef std::vector<Lit> CL;
  pepin.add_clause(CL{itol(1), itol(2), itol(3)});
  pepin.add_clause(CL{itol(1), itol(-5)});
  auto weigh_num_sols = pepin.get_low_prec_appx_weighted_sol();
  std::cout << "Solution: " << std::scientific << std::setprecision(10)
    << *weigh_num_sols << std::endl;

  return 0;
}

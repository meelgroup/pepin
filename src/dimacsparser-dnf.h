/*****************************************************************************
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
CryptoMiniSat -- Copyright (C) 2009-2020 Authors of CryptoMiniSat, see AUTHORS file
Pepin -- Copyright (C) 2021 -- Mate Soos, Kuldeep S. Meel

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
******************************************************************************/

#ifndef DIMACSPARSER_H
#define DIMACSPARSER_H

#include <string.h>
#include "streambuffer.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <cassert>
#include "pepin.h"

using std::vector;
using std::cout;
using std::endl;
using PepinNS::Lit;

template <class C, class S>
class DimacsParser
{
    public:
        DimacsParser(S* solver, unsigned _verbosity);

        template <class T> bool parse_DIMACS(
            T input_stream);
        uint64_t max_var = std::numeric_limits<uint64_t>::max();
        vector<uint32_t> sampling_vars;
        vector<double> weights;
        const std::string dimacs_spec = "http://www.satcompetition.org/2009/format-benchmarks2009.html";
        const std::string please_read_dimacs = "\nPlease read DIMACS specification at http://www.satcompetition.org/2009/format-benchmarks2009.html";

    private:
        bool parse_DIMACS_main(C& in);
        bool readClause(C& in);
        bool parse_and_add_clause(C& in);
        bool match(C& in, const char* str);
        bool parse_header(C& in);
        std::string stringify(uint32_t x) const;
        bool parseWeight(C& in);

        S* solver;
        unsigned verbosity;

        //Stat
        size_t lineNum;

        //check header strictly
        bool header_found = false;
        int num_header_vars = 0;
        int num_header_cls = 0;

        //Reduce temp overhead
        vector<Lit> lits;

        size_t norm_clauses_added = 0;
};

template<class C, class S>
DimacsParser<C, S>::DimacsParser(
    S* _solver
    , unsigned _verbosity
):
    solver(_solver)
    , verbosity(_verbosity)
    , lineNum(0)
{
}

template<class C, class S>
std::string DimacsParser<C, S>::stringify(uint32_t x) const
{
    std::ostringstream o;
    o << x;
    return o.str();
}

template<class C, class S>
bool DimacsParser<C, S>::readClause(C& in)
{
    int32_t parsed_lit;
    uint32_t var;
    for (;;) {
        if (!in.parseInt(parsed_lit, lineNum)) {
            return false;
        }
        if (parsed_lit == 0) {
            break;
        }

        var = std::abs(parsed_lit)-1;

        if (var > max_var) {
            std::cerr
            << "ERROR! "
            << "Variable requested is too large for DIMACS parser parameter: "
            << var << endl
            << "--> At line " << lineNum+1
            << please_read_dimacs
            << endl;
            return false;
        }

        if (var >= (1ULL<<28)) {
            std::cerr
            << "ERROR! "
            << "Variable requested is far too large: " << var + 1 << endl
            << "--> At line " << lineNum+1
            << please_read_dimacs
            << endl;
            return false;
        }

        if (!header_found) {
            std::cerr
            << "ERROR! "
            << "DIMACS header ('p cnf vars cls') never found!" << endl;
            return false;
        }

        if ((int)var >= num_header_vars) {
            std::cerr
            << "ERROR! "
            << "Variable requested is larger than the header told us." << endl
            << " -> var is : " << var + 1 << endl
            << " -> header told us maximum will be : " << num_header_vars << endl
            << " -> At line " << lineNum+1
            << endl;
            return false;
        }

        lits.push_back( (parsed_lit > 0) ? Lit(var, false) : Lit(var, true) );
        if (*in != ' ') {
            std::cerr
            << "ERROR! "
            << "After last element on the line must be 0" << endl
            << "--> At line " << lineNum+1
            << please_read_dimacs
            << endl
            << endl;
            return false;
        }
    }

    return true;
}

template<class C, class S>
bool DimacsParser<C, S>::match(C& in, const char* str)
{
    for (; *str != 0; ++str, ++in)
        if (*str != *in)
            return false;
    return true;
}

template<class C, class S>
bool DimacsParser<C, S>::parseWeight(C& in)
{
    if (match(in, "w ")) {
        int32_t slit;
        uint32_t dividend;
        uint32_t divisor;
        if (in.parseInt(slit, lineNum)
            && in.parseRatio(dividend, divisor, lineNum)
        ) {
            if (slit == 0) {
                cout << "ERROR: Cannot define weight of literal 0!" << endl;
                exit(-1);
            }
            uint32_t var = std::abs(slit)-1;
            solver->set_var_weight(var, dividend, divisor);
            return true;
        } else {
            cout << "ERROR: weight is incorrect on line " << lineNum << endl;
            exit(-1);
        }
    } else {
        cout << "ERROR: weight is not given on line " << lineNum << endl;
        exit(-1);
    }
    return true;
}

template<class C, class S>
bool DimacsParser<C, S>::parse_header(C& in)
{
    ++in;
    in.skipWhitespace();
    std::string str;
    in.parseString(str);
    if (str == "cnf") {
        cout << "WARN: header says it's a CNF not a DNF!" << endl;
        //exit(-1);
    }
    if (str == "dnf" || str == "cnf") {
        if (header_found) {
            std::cerr << "ERROR: CNF header ('p dnf vars cls') found twice in file! Exiting." << endl;
            exit(-1);
        }
        header_found = true;

        if (!in.parseInt(num_header_vars, lineNum)
            || !in.parseInt(num_header_cls, lineNum)
        ) {
            return false;
        }
        if (verbosity) {
            cout << "c -- header says num vars:   " << std::setw(12) << num_header_vars << endl;
            cout << "c -- header says num clauses:" <<  std::setw(12) << num_header_cls << endl;
        }
        if (num_header_vars < 0) {
            std::cerr << "ERROR: Number of variables in header cannot be less than 0" << endl;
            return false;
        }
        if (num_header_cls < 0) {
            std::cerr << "ERROR: Number of clauses in header cannot be less than 0" << endl;
            return false;
        }
        solver->set_n_cls(num_header_cls);
        solver->new_vars(num_header_vars);
    } else {
        std::cerr
        << "PARSE ERROR! Unexpected char (hex: " << std::hex
        << std::setw(2)
        << std::setfill('0')
        << "0x" << *in
        << std::setfill(' ')
        << std::dec
        << ")"
        << " At line " << lineNum+1
        << "' in the header"
        << please_read_dimacs
        << endl;
        return false;
    }

    return true;
}

template<class C, class S>
bool DimacsParser<C, S>::parse_and_add_clause(C& in)
{
    lits.clear();
    if (!readClause(in)) {
        return false;
    }
    in.skipWhitespace();
    if (!in.skipEOL(lineNum)) {
        return false;
    }
    lineNum++;
    solver->add_clause(lits, norm_clauses_added);
    norm_clauses_added++;
    return true;
}

template<class C, class S>
bool DimacsParser<C, S>::parse_DIMACS_main(C& in)
{
    std::string str;

    for (;;) {
        in.skipWhitespace();
        switch (*in) {
        case EOF:
            return true;
        case 'p':
            if (!parse_header(in)) {
                return false;
            }
            in.skipLine();
            lineNum++;
            break;
        case 'w':
            if (!parseWeight(in)) {
                return false;
            }
            in.skipLine();
            lineNum++;
            break;
        case 'c':
            ++in;
            in.skipLine();
            lineNum++;
            break;
        case '\n':
            if (verbosity) {
                std::cout
                << "c WARNING: Empty line at line number " << lineNum+1
                << " -- this is not part of the DIMACS specifications ("
                << dimacs_spec << "). Ignoring."
                << endl;
            }
            in.skipLine();
            lineNum++;
            break;
        default:
            if (!parse_and_add_clause(in)) {
                return false;
            }
            break;
        }
    }

    return true;
}

template <class C, class S>
template <class T>
bool DimacsParser<C, S>::parse_DIMACS(
    T input_stream)
{
    assert(solver->nVars() == 0);

    C in(input_stream);
    if ( !parse_DIMACS_main(in)) {
        return false;
    }

    if (verbosity) {
        cout
        << "c -- clauses added: " << norm_clauses_added << endl;
    }

    return true;
}

#endif //DIMACSPARSER_H

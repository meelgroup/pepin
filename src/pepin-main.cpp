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

#include <iostream>
#include <iomanip>
#include <string>
#include "argparse.hpp"

#include "time_mem.h"
#include "pepin.h"
#include "dimacsparser-dnf.h"

// To allow breaking on division by zero etc
#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;

double startTime;
int verb = 1;
int seed = 1;
double epsilon = 0.5;
double delta = 0.36;
PepinNS::Pepin* dnfs;
int force_eager = false;
int fast_center_calc = 1;

#define myopt(name, var, fun, hhelp) \
    program.add_argument(name) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)
#define myopt2(name1, name2, var, fun, hhelp) \
    program.add_argument(name1, name2) \
        .action([&](const auto& a) {var = std::fun(a.c_str());}) \
        .default_value(var) \
        .help(hhelp)

argparse::ArgumentParser program = argparse::ArgumentParser("arjun", PepinNS::Pepin::get_version_sha1(),
        argparse::default_arguments::help);

void print_version() {
    cout << "c [dnfs] SHA revision: " << PepinNS::Pepin::get_version_sha1() << endl;
    cout << "c [dnfs] Compilation environment: " << PepinNS::Pepin::get_compilation_env() << endl;
    std::exit(0);
}

void add_dnfs_options() {
    myopt2("-v", "--verb", verb, atoi, "Verbosity");
    program.add_argument("-v", "--version") \
        .action([&](const auto&) {print_version(); exit(0);}) \
        .flag()
        .help("Print version and exit");
    myopt2("--epsilon","-e", epsilon, atof, "epsilon");
    myopt2("--delta","-d", delta, atof, "delta");
    myopt2("--verb","-v", verb, atoi, "verbosity");
    myopt2("--seed","-s", seed, atoi, "Seed");
    myopt("--eager", force_eager, atoi, "Force eager");
    myopt("--fastcenter", fast_center_calc, atoi, "fast center calculation");

    program.add_argument("files").remaining().help("input file and output file");
}

void readInAFile(const string& filename)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN>, PepinNS::Pepin> parser(dnfs, verb);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ>, PepinNS::Pepin> parser(dnfs, verb);
    #endif

    if (in == NULL) {
        std::cerr
        << "ERROR! Could not open file '"
        << filename
        << "' for reading: " << strerror(errno) << endl;

        std::exit(-1);
    }

    if (!parser.parse_DIMACS(in)) exit(-1);

    #ifndef USE_ZLIB
        fclose(in);
    #else
        gzclose(in);
    #endif
}

int main(int argc, char** argv)
{
    // Die on division by zero etc.
    #if defined(__GNUC__) && defined(__linux__)
    feenableexcept(FE_INVALID   |
                   FE_DIVBYZERO |
                   FE_OVERFLOW
                  );
    #endif

    //Reconstruct the command line so we can emit it later if needed
    string command_line;
    for(int i = 0; i < argc; i++) {
        command_line += string(argv[i]);
        if (i+1 < argc) {
            command_line += " ";
        }
    }
    add_dnfs_options();
    try {
        program.parse_args(argc, argv);
        if (program.is_used("--help")) {
            cout
            << "DNF probabilistic aproximate counter." << endl << endl
            << "pepin inputfile" << endl;
            cout << program << endl;
            std::exit(0);
        }
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        exit(-1);
    }
    dnfs = new PepinNS::Pepin(epsilon, delta, seed, verb);
    dnfs->set_force_eager(force_eager);
    dnfs->set_fast_center_calc(fast_center_calc);

    cout << "c [dnfs] Pepin Version: " << dnfs->get_version_sha1() << endl;
    if (verb >= 2) {
        cout << "c [dnfs] compilation environment: " << dnfs->get_compilation_env()
        << endl;

        cout
        << "c [dnfs] executed with command line: "
        << command_line
        << endl;
    }

    double starTime = cpuTime();
    cout << "c [dnfs] using epsilon: " << epsilon
        << " delta: "<< delta << " seed: " << seed << endl;

    //parsing the input
    vector<std::string> files;
    try {
        files = program.get<std::vector<std::string>>("files");
        if (files.size() != 1) {
            cout << "ERROR: you must pass at exactly file: an INPUT file" << endl;
            exit(-1);
        }
    } catch (std::logic_error& e) {
        cout << "ERROR: you must give an input file" << endl;
        exit(-1);
    }
    readInAFile(files[0]);

    auto low_prec_num_points = dnfs->get_low_prec_appx_num_points();
    cout << "c [dnfs] Low-precision approx num points: " << std::fixed << std::setprecision(0)
            << *low_prec_num_points << std::setprecision(10) << endl;
    auto weigh_num_sols = dnfs->get_appx_weighted_sol();
    cout << "c [dnfs] Weight no. solutions: " << *weigh_num_sols << endl;

    auto low_prec_weigh_num_sols = dnfs->get_low_prec_appx_weighted_sol();
    cout << "c [dnfs] Low-precision weighted no. solutions: " << std::scientific << std::setprecision(4) << *low_prec_weigh_num_sols << endl;

    cout << "c [dnfs] finished T: " << std::setprecision(2) << std::fixed << (cpuTime() - starTime)
    << endl;

    delete dnfs;
    return 0;
}

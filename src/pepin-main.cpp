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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iomanip>
#include <vector>
#include <atomic>
#include <fstream>
#include <sstream>
#include <string>
#include <signal.h>

#include "time_mem.h"
#include "dimacsparser-dnf.h"
#include "pepin.h"

// To allow breaking on division by zero etc
#if defined(__GNUC__) && defined(__linux__)
#include <fenv.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;

po::options_description dnfs_options = po::options_description("Pepin options");
po::options_description help_options;
po::variables_map vm;
po::positional_options_description p;
double startTime;
int verb = 0;
int seed = 1;
double epsilon = 0.5;
double delta = 0.36;
Pepin* dnfs;
int force_eager = false;
int fast = 1;

// static void signal_handler(int) {
//     cout << endl << "c [dnfs] INTERRUPTING ***" << endl << std::flush;
//     common.interrupt_asap = true;
// }

void add_dnfs_options()
{
    dnfs_options.add_options()
    ("help,h", "Prints help")
    ("version", "Print version info")
    ("input", po::value<string>(), "file to read")
    ("verb,v", po::value(&verb)->default_value(verb), "verbosity")
    ("seed,s", po::value(&seed)->default_value(seed), "Seed")
    ("epsilon,e", po::value(&epsilon)->default_value(epsilon), "epsilon")
    ("eager", po::value(&force_eager)->default_value(force_eager), "Force eager")
    ("delta,d", po::value(&delta)->default_value(delta), "delta")
    ("fast", po::value(&fast)->default_value(fast), "fast")
    ;

    help_options.add(dnfs_options);
}

void add_supported_options(int argc, char** argv)
{
    add_dnfs_options();
    p.add("input", 1);

    try {
        po::store(po::command_line_parser(argc, argv).options(help_options).positional(p).run(), vm);
        if (vm.count("help"))
        {
            cout
            << "Minimal projection set finder" << endl;

            cout
            << "dnfs [options] inputfile" << endl << endl;

            cout << help_options << endl;
            std::exit(0);
        }

        if (vm.count("version")) {
            dnfs = new Pepin(epsilon, delta, seed);
            cout << "c [dnfs] SHA revision: " << dnfs->get_version_info() << endl;
            cout << "c [dnfs] Compilation environment: " << dnfs->get_compilation_env() << endl;
            delete dnfs;
            std::exit(0);
        }

        po::notify(vm);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::unknown_option> >& c
    ) {
        cerr
        << "ERROR: Some option you gave was wrong. Please give '--help' to get help" << endl
        << "       Unkown option: " << c.what() << endl;
        std::exit(-1);
    } catch (boost::bad_any_cast &e) {
        std::cerr
        << "ERROR! You probably gave a wrong argument type" << endl
        << "       Bad cast: " << e.what()
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_option_value> >& what
    ) {
        cerr
        << "ERROR: Invalid value '" << what.what() << "'" << endl
        << "       given to option '" << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::multiple_occurrences> >& what
    ) {
        cerr
        << "ERROR: " << what.what() << " of option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::required_option> >& what
    ) {
        cerr
        << "ERROR: You forgot to give a required option '"
        << what.get_option_name() << "'"
        << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::too_many_positional_options_error> >& what
    ) {
        cerr
        << "ERROR: You gave too many positional arguments. Only the input CNF can be given as a positional option." << endl;
        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::ambiguous_option> >& what
    ) {
        cerr
        << "ERROR: The option you gave was not fully written and matches" << endl
        << "       more than one option. Please give the full option name." << endl
        << "       The option you gave: '" << what.get_option_name() << "'" <<endl
        << "       The alternatives are: ";
        for(size_t i = 0; i < what.alternatives().size(); i++) {
            cout << what.alternatives()[i];
            if (i+1 < what.alternatives().size()) {
                cout << ", ";
            }
        }
        cout << endl;

        std::exit(-1);
    } catch (boost::exception_detail::clone_impl<
        boost::exception_detail::error_info_injector<po::invalid_command_line_syntax> >& what
    ) {
        cerr
        << "ERROR: The option you gave is missing the argument or the" << endl
        << "       argument is given with space between the equal sign." << endl
        << "       detailed error message: " << what.what() << endl
        ;
        std::exit(-1);
    }
}

void readInAFile(const string& filename)
{
    #ifndef USE_ZLIB
    FILE * in = fopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<FILE*, FN>, Pepin> parser(dnfs, verb);
    #else
    gzFile in = gzopen(filename.c_str(), "rb");
    DimacsParser<StreamBuffer<gzFile, GZ>, Pepin> parser(dnfs, verb);
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
    add_supported_options(argc, argv);
    dnfs = new Pepin(epsilon, delta, seed, verb);
    dnfs->set_force_eager(force_eager);
    dnfs->set_fast(fast);

    cout << "c Pepin Version: " << dnfs->get_version_info() << endl;
    cout << "c compilation environment: " << dnfs->get_compilation_env()
    << endl;

    cout
    << "c executed with command line: "
    << command_line
    << endl;

    double starTime = cpuTime();
    cout << "c [dnfs] using seed: " << seed << endl;

    //parsing the input
    if (vm.count("input") == 0) {
        cout << "ERROR: you must pass a file" << endl;
        exit(-1);
    }
    const string inp = vm["input"].as<string>();
    readInAFile(inp);

    cout << "c [dnfs] finished T: " << std::setprecision(2) << std::fixed << (cpuTime() - starTime)
    << endl;

    delete dnfs;
    return 0;
}

#!/usr/bin/python3

#  Pepin
#
#  Copyright (c) 2021, Mate Soos and Kuldeep S. Meel. All rights reserved.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.

import os
import subprocess
import decimal
import time
import random
import optparse
import glob

def shuffle(fname):
    tmpfname = "tmp.dnf"
    with open(fname, "r") as f1:
        with open(tmpfname, "w") as f2:
            buf = []
            for line in f1:
                if "p cnf" in line or "p dnf" in line:
                    f2.write(line)
                    continue
                line = line.strip()
                if len(line) == 0:
                    continue
                line = line.split()
                assert len(line) > 0
                line = line[:len(line)-1]
                random.shuffle(line)
                buf.append(line)

            random.shuffle(buf)

            for b in buf:
                f2.write(" ".join(b)+" 0\n")
    return tmpfname


def get_rel_count_DNFKLM(fullfilename, nvars, seed_here, options):
    out = subprocess.check_output(["./DNFKLM", options.epsilon, options.delta,
                                   fullfilename, "%s" % seed_here])
    count = None
    for l in out.splitlines():
        l = l.decode()
        if "Final count is" in l:
            #print(l)
            count = l.split(" ")[3]
            count_left = count.split("^")[0]
            assert count_left == "2"
            count_right = count.split("^")[1]
            #print(count_left, " -- ", count_right)
            count = decimal.Decimal(2) ** decimal.Decimal(count_right)
            break;
    if count is None:
        print("Error, DNFKLM didn't return a solution")
        exit(-1)

    #print("DNFKLM count:", count)

    relcount = count/(decimal.Decimal(2)**decimal.Decimal(nvars))
    #print("DNFKLM Relative count:" , relcount)

    return relcount


def get_rel_count_Pepin(fullfilename, seed_here, options, sparse=None):
    if sparse is None:
        sparse = options.sparse
    run_command = ["./pepin",
                   "--epsilon", options.epsilon, "--delta", options.delta,
                   "--seed", "%s" % seed_here,
                   "--sparse", "1" if sparse else "0",
                   fullfilename]
    if options.verbose:
        print("Running Pepin (sparse=%s) on file %s" % (sparse, fullfilename))
        print("Running: %s" % " ".join(run_command))
    out2 = subprocess.check_output(run_command)

    relcount2 = None
    for l in out2.splitlines():
        l = l.decode()
        if "Low-precision weighted no. solutions" in l:
            #print(l)
            relcount2 = decimal.Decimal(l.split()[6])

    if relcount2 is None:
        print("Error, Pepin failed to count!")
        exit(-1)

    #print("Pepin Relative count:" , relcount2)
    return relcount2


def get_num_vars(fullfilename):
    nvars = None
    with open(fullfilename, "r") as f:
        for line in f:
            line = line.strip()
            if "p dnf" in line:
                nvars = line.split(" ")[2]

    assert nvars is not None
    return nvars


def print_decimal(x):
    return x
    #return x.normalize()


def print_quant_perc(x):
    return decimal.Decimal(100.0)-(x).quantize(decimal.Decimal('.001'))


def test(fname, num_tested, options):
    fullfilename = "%s/%s" % (options.dir, fname)
    nvars = get_num_vars(fullfilename)
    if options.verbose:
        print("----------------------------------------- file num: %-10d" % num_tested)
        print("Dealing with file %s/%s" % (options.dir, fname))
        print("Num vars:", nvars)

    if options.self_test:
        # checking smallest vs largest DNFKLM
        t = time.time()
        relcounts = []
        for i in range(options.self_count):
            newf = shuffle(fullfilename)
            relcounts.append(get_rel_count_DNFKLM(newf, nvars, i, options))
        s = sorted(relcounts)

        if options.verbose:
            print("Smallest vs Largest DNFKLM over %d runs: %s %%"
                   % (options.self_count, print_quant_perc(print_decimal(s[0]/s[len(s)-1])*100)),
                  " T %-3.3f" % (time.time()-t))

        if s[0]/s[len(s)-1] < decimal.Decimal(0.5):
            print("WARNING: too large difference by the same counter, DNFKLM!!")
            #exit(-1)

        # checking smallest vs largest DNF Stream
        t = time.time()
        relcounts = []
        for i in range(options.self_count):
            newf = shuffle(fullfilename)
            relcounts.append(get_rel_count_Pepin(newf, i, options))
        s = sorted(relcounts)

        if options.verbose:
            print("Smallest vs Largest count by Pepin over  %d runs: %s %%"
                  % (options.self_count, print_quant_perc(print_decimal(s[0]/s[len(s)-1])*100)),
                  " T %-3.3f" % (time.time()-t))
        if s[0]/s[len(s)-1] < decimal.Decimal(0.5):
            print("WARNING: too large difference by the same counter, Pepin!!")
            #exit(-1)

    # compare DNFKLM against Pepin in both dense and sparse modes
    relcount = get_rel_count_DNFKLM(fullfilename, nvars, options.seed, options)
    for sparse_mode in (False, True):
        relcount2 = get_rel_count_Pepin(fullfilename, options.seed, options, sparse=sparse_mode)
        mode_label = "sparse" if sparse_mode else "dense"
        if options.verbose:
            print("Relcount DNFKLM:           ", relcount)
            print("Relcount Pepin (%s):  " % mode_label, relcount2)
        print("Diff DNFKLM vs Pepin(%s):  %s %%"
              % (mode_label,
                 print_quant_perc(print_decimal(relcount/relcount2)*100)))
        if relcount/relcount2 > 3 or relcount2/relcount > 3:
            print("!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!")
            print("Something is WRONG: DNFKLM vs Pepin (%s)" % mode_label)
            print("Diff:", print_decimal(relcount/relcount2))
            print("DNFKLM count:" , print_decimal(relcount))
            print("DNF Streaming count (%s):" % mode_label, print_decimal(relcount2))
            exit(-1)


def run_tests():
    print("Running tests in dir %s..." % options.dir)
    files = os.listdir(options.dir)
    num_tested = 0
    for f in files:
        if ".dnf" in f:
            num_tested+=1
            test(f, num_tested, options)
    for file in glob.glob("tests/*.dnf"):
        os.unlink(file)
    print("Finished all OK! Tested %d Pepin" % num_tested)


class PlainHelpFormatter(optparse.IndentedHelpFormatter):
    def format_description(self, description):
        if description:
            return description + "\n"
        else:
            return ""

def set_up_parser():
    parser = optparse.OptionParser(usage="Run with default options",
                                   description="usage: %prog [options]",
                                   formatter=PlainHelpFormatter())

    parser.add_option("--verbose", "-v", action="store_true", default=False,
                      dest="verbose", help="Print more output")
    parser.add_option("--seed", dest="seed", default=1,
                      help="Fuzz test start seed. Otherwise, random seed is picked"
                      " (printed to console)", type=int)
    parser.add_option("--self", dest="self_test", default=False,
                      help="Perform self-test", action="store_true")
    parser.add_option("--selfcount", dest="self_count", default=100,
                      type=int, help="How many tests with self-test")
    parser.add_option("--dir", dest="dir", type=str,
                      default="tests",
                      help="Tests are in this dir")
    parser.add_option("--instances", dest="instances", default="3",
                      type=str, help="How many files per type")
    parser.add_option("--epsilon", dest="epsilon", default="0.5",
                      type=str, help="Epsilon to use")
    parser.add_option("--delta", dest="delta", default="0.35",
                      type=str, help="Delta to use")
    parser.add_option("--sparse", dest="sparse", default=False,
                      action="store_true", help="Use sparse representation (default: dense)")
    return parser


if __name__ == "__main__":
    # setup
    decimal.getcontext().prec = 300

    # options
    parser = set_up_parser()
    (options, args) = parser.parse_args()
    random.seed(options.seed)

    # run tests
    os.system("rm -f tests/rand*.dnf")
    os.makedirs(options.dir, exist_ok=True)

    t = time.time();
    print("Generating small random instances...")
    ret = os.system("./random_dnf_generator.py --instances %s -n '101,1002,301' --cldensity '0.2,1.1,0.2' --clsize '10,80,20' --noscaling --monotone %s/" % (options.instances, options.dir))
    assert ret == 0, "random DNF generation failed"
    print("Done, T: ", time.time()-t)
    files = os.listdir(options.dir)
    run_tests()

    print("Generating large random instances...")
    t = time.time();
    ret = os.system("./random_dnf_generator.py --instances %s -n '1000,10001,3000' --cldensity '0.01,0.04,0.01' --clsize '100,350,50' --noscaling --monotone %s/" % (options.instances, options.dir))
    assert ret == 0, "random DNF generation failed"
    print("Done, T: ", time.time()-t)
    run_tests()

    print("Generating random instances with odd num vars...")
    t = time.time();
    ret = os.system("./random_dnf_generator.py --instances %s -n '100,110,1' --cldensity '0.1,0.5,0.2' --clsize '10,100,20' --noscaling --monotone %s/" % (options.instances, options.dir))
    assert ret == 0, "random DNF generation failed"
    print("Done, T: ", time.time()-t)
    run_tests()

    print("Generating random instances with VERY large num vars...")
    t = time.time();
    ret = os.system("./random_dnf_generator.py --instances %s -n 200000,200001,200000 --cldensity 0.1,0.2,0.1 --clsize 10,20,50 --noscaling --monotone %s/" % (options.instances, options.dir))
    assert ret == 0, "random DNF generation failed"
    print("Done, T: ", time.time()-t)
    run_tests()

    # Many-clauses regime — the sparse target. The scaffolded generator
    # is O(N^2) on this shape, so use fast_gen.py instead.
    print("Generating MANY-clauses instances (sparse target)...")
    os.system("rm -f %s/rand*.dnf %s/big_*.dnf" % (options.dir, options.dir))
    t = time.time()
    instances = int(options.instances)
    for idx in range(instances):
        # Each instance: 10k vars x 50k clauses x 3-lit. Big enough to
        # expose dense vs sparse perf gap, small enough for DNFKLM.
        out = "%s/big_%d.dnf" % (options.dir, idx)
        ret = os.system("./fast_gen.py 10000 50000 3 %s %d" % (out, idx + 1))
        assert ret == 0, "fast_gen.py failed"
    print("Done, T: ", time.time()-t)
    run_tests()

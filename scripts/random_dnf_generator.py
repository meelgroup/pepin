#!/usr/bin/env python3

#  Copyright (c) 2019 Kuldeep S. Meel, Aditya A. Shrotri, Moshe Y. Vardi
#                2021 Updates by Mate Soos
#                All rights reserved.
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


import random
from optparse import OptionParser
random.seed(1)

usage = "usage: %prog [options] outdir"
parser = OptionParser(usage)
parser.add_option("--instances", dest="instances", type=int,
                  help="How many instances to generate")
parser.add_option("--n", dest="n_info", default="",
                  help="start-stop-step for N (number of variables) in the form of 'start,stop,end'")
parser.add_option("--mdensity", dest="mdensity_info", default="",
                  help="start-stop-step for M density (i.e. number of clauses) in the form of 'start,stop,end'")
parser.add_option("--msize", dest="msize_info", default="",
                  help="start-stop-step for clause size in the form of 'start,stop,end'")
parser.add_option("--scaling",
                  action="store_true", dest="scaling", default=False,
                  help="Scaling size with N")
parser.add_option("--noscaling",
                  action="store_true", dest="noscaling", default=False,
                  help="Do not scale size with N")
parser.add_option("--monotone",
                  action="store_true", dest="monotone", default=False,
                  help="Monotone")
parser.add_option("--nomonotone",
                  action="store_true", dest="nomonotone", default=False,
                  help="Not monotone")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="Print more to stdout")

(options, args) = parser.parse_args()
if len(args) != 1:
    parser.error("ERROR: You must pass EXACTLY one argument, the output dir")
    exit(-1)
outputDir = args[0]

if options.instances is None:
    print("ERROR: you must give the number of instances to generate");
    exit(-1)

# parse scaling
if options.scaling and options.noscaling:
    print("ERROR: Cannot both have scaling and non-scaling options")
    exit(-1)

if (options.scaling or options.noscaling) == False:
    print("ERROR: Must select either --scaling or --noscaling")
    exit(-1)

mSizeType = 1
if options.scaling:
    mSizeType = 0

#parse monotone
if options.monotone and options.nomonotone:
    print("ERROR: Cannot both have monotone and non-monotone options")
    exit(-1)

if (options.monotone or options.nomonotone) == False:
    print("ERROR: Must select either --monotone or --nomonotone")
    exit(-1)

monotone = 0
if options.monotone:
    monotone = 1


def parse_triplet(input :str, name :str, func):
    if input.strip() == "":
        print("ERROR: you must give %s option in the form of 'start, stop, step'" % name);
        exit(-1)

    input2 = input.strip().split(',')
    if len(input2) != 3:
        print("Error: cannot parse %s option, it doesn't have 3 parts. You must give 'start,end,step'" % name);
        exit(-1);

    start : int = 0
    stop :int  = 0
    step :int = 0
    try:
        start = func(input2[0])
    except:
        print("ERROR: Start value of %s is not an integer" % name)
        exit(-1)
    try:
        stop = func(input2[1])
    except:
        print("ERROR: Stop value of %s is not an integer" % name)
        exit(-1)

    try:
        step = func(input2[2])
    except:
        print("ERROR: Step value of %s is not an integer" % name)
        exit(-1)

    if start > stop:
        print("ERROR: 'start' is larger than 'stop' for option %s" % name)
        exit(-1)

    return start, stop, step

nLow, nHigh, nStep = parse_triplet(options.n_info, "--n", int)
mDensityLow, mDensityHigh, mDensityStep = parse_triplet(options.mdensity_info, "--mdensity", float)
mSizeLow, mSizeHigh, mSizeStep  = parse_triplet(options.msize_info, "--msize", float)

for n in range(nLow,nHigh,nStep):
    m = mDensityLow
    while m < mDensityHigh:
        mSizeLow1 = mSizeLow
        mSizeHigh1 = mSizeHigh
        mSizeStep1 = mSizeStep
        if mSizeType==0:
            mSizeLow1 = n*mSizeLow
            mSizeHigh1 = n*mSizeHigh
            mSizeStep1 = n*mSizeStep
            print(n,mSizeLow,mSizeHigh,mSizeStep)
        for k in range(int(mSizeLow1),int(mSizeHigh1),int(mSizeStep1)):
            for i in range(options.instances):
                opStr = "p dnf "+str(n)+" "+str(int(n*m))+'\n'
                for j in range(int(n*m)):
                    sampVars = random.sample(range(1,n+1),k)
                    if (j%200 == 0):
                        sampVars = random.sample(range(1,n+1),min(n, 300*k))
                    for l in sampVars:
                        if monotone == 1:
                            opStr += str(l)+" "
                        else:
                            if random.random()>0.5:
                                opStr += str(l)+" "
                            else:
                                opStr += "-"+str(l)+" "
                    opStr += "0\n"
                f = open(outputDir+"/randomDNF_"+str(n)+"_"+str(int(n*m))+"_"+str(k)+"_"+str(i)+".dnf",'w')
                f.write(opStr)
                f.close()

        m += mDensityStep

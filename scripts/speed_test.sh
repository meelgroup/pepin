#!/bin/bash
set -x
rm -f tests/*.dnf
./random_dnf_generator.py --instances 1 -n '200000,200001,200000' --cldensity '0.1,0.2,0.1' --clsize '10,20,50' --noscaling --monotone tests
FILE="tests/randomDNF_200000_20000_10_0.dnf"
./pepin $FILE | tail -n 3
./pepin --sparse 1 $FILE | tail -n 3

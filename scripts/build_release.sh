#!/bin/bash

set -e

rm -rf cm* CM* lib* Testing* tests* include tests *.cmake Make* dnfstream src countbysample hypervolstream*
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j26
# make test

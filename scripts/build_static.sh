#!/bin/bash

set -e

rm -rf cm* CM* lib* Testing* tests* *.cmake Make* dnf_streaming_count countbysample src
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DSTATICCOMPILE=ON ..
make -j26

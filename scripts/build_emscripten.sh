#!/bin/bash

set -e
rm -rf cm* CM* lib* cryptomini* Testing* tests* pycryptosat include tests utils Make*
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DENABLE_TESTING=OFF ..
emmake make -j26
emmake make install
cp pepin.wasm ../html
cp $EMINSTALL/bin/pepin.js ../html

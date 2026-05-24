#!/bin/bash
set -euxo pipefail

rsync -vaP pepin.js pepin.wasm  index.html msoos.org:/var/www/pepin/

#!/bin/bash
set -euxo pipefail

make clean
make -j main

mkdir -p out

(time ./main bigl w   2147483647  103 2) &> out/l-100.txt
(time ./main bigl w   2147483647  211 2) &> out/l-200.txt
(time ./main bigl w   2147483647  419 2) &> out/l-400.txt
(time ./main bigl w   2147483647  811 2) &> out/l-800.txt
(time ./main bigl w   2147483647 1019 2) &> out/l-1000.txt
(time ./main bigc w   2147483647 1019 2) &> out/c-1000.txt
(time ./main bigc w   2147483647 2003 2) &> out/c-2000.txt
(time ./main bigc w   2147483647 4003 2) &> out/c-4000.txt
(time ./main bigc w   2147483647 8011 2) &> out/c-8000.txt


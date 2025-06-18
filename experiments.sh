#!/bin/bash
set -euxo pipefail

make clean
make -j main

mkdir -p out

(time ./main bigl w   2147483647  211 2) &> out/l-200.txt
(time ./main bigc w   2147483647  211 2) &> out/c-200.txt
(time ./main bigl w   2147483647  419 2) &> out/l-400.txt
(time ./main bigc w   2147483647  419 2) &> out/c-400.txt
(time ./main bigl w   2147483647  607 2) &> out/l-600.txt
(time ./main bigc w   2147483647  607 2) &> out/c-600.txt
(time ./main bigl w   2147483647  811 2) &> out/l-800.txt
(time ./main bigc w   2147483647  811 2) &> out/c-800.txt
(time ./main bigl w   2147483647 1019 2) &> out/l-1000.txt
(time ./main bigc w   2147483647 1019 2) &> out/c-1000.txt
(time ./main bigl w   2147483647 2003 2) &> out/l-2000.txt
(time ./main bigc w   2147483647 2003 2) &> out/c-2000.txt
(time ./main bigl w   2147483647 3011 2) &> out/l-3000.txt
(time ./main bigc w   2147483647 3011 2) &> out/c-3000.txt
(time ./main bigl w   2147483647 4003 2) &> out/l-4000.txt
(time ./main bigc w   2147483647 4003 2) &> out/c-4000.txt
(time ./main bigl w   2147483647 5003 2) &> out/l-5000.txt
(time ./main bigc w   2147483647 5003 2) &> out/c-5000.txt
(time ./main bigl w   2147483647 6007 2) &> out/l-6000.txt
(time ./main bigc w   2147483647 6007 2) &> out/c-6000.txt
(time ./main bigl w   2147483647 7019 2) &> out/l-7000.txt
(time ./main bigc w   2147483647 7019 2) &> out/c-7000.txt
(time ./main bigl w   2147483647 8011 2) &> out/l-8000.txt
(time ./main bigc w   2147483647 8011 2) &> out/c-8000.txt
(time ./main bigl w   2147483647 9007 2) &> out/l-9000.txt
(time ./main bigc w   2147483647 9007 2) &> out/c-9000.txt
(time ./main bigl w   2147483647 10007 2) &> out/l-10000.txt
(time ./main bigc w   2147483647 10007 2) &> out/c-10000.txt
(time ./main bigl w   2147483647 11681 2) &> out/l-11681.txt
(time ./main bigc w   2147483647 11681 2) &> out/c-11681.txt


#!/usr/bin/env sh
make
dist/osx/debug/test.app/Contents/MacOS/test discreteLaplacian > out.txt
gnuplot -p -e "set term wx; set style data lines; splot 'out.txt' using 1:2:3; pause -1"

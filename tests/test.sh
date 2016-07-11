#!/usr/bin/env sh
make clean all
dist/osx/debug/test.app/Contents/MacOS/test discreteLaplacian > out.txt
gnuplot -e "set style data lines; splot 'out.txt' using 1:2:3"

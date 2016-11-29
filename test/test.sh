#!/usr/bin/env sh
#dist/osx/debug/test.app/Contents/MacOS/test discreteLaplacian > out.txt
lua -lmake && dist/linux/debug/test discreteLaplacian > out.txt && ./plot_results.gnuplot

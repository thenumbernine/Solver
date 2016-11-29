#! /usr/bin/gnuplot -p

set style data lines
set terminal png size 800,600

set log y
set output 'solver_residual.png'
plot 'solver.txt' using 1:2 title 'residual'
unset log y

set output 'phi.png'
splot 'out.txt' using 1:2:3 title 'phi'

# might error if i'm not running jfnk ...
set log y
set output 'gmres_convergence.png'
plot 'gmres.txt' using 2:3:1 palette
unset log y

# might error if solver.txt has no alpha ... (ie the linear solvers)
set output 'solver_alpha.png'
plot 'solver.txt' using 1:3 title 'alpha'

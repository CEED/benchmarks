#! /usr/bin/gnuplot

set xlabel "DOFS" font ",14"
set ylabel "DOFS/s" font ",14"
set logscale x 2
set logscale y 3
set term png

set title  "Nek5000 - BP1 - Vector" font ",16"
plot for [i=start:end] "lx".i."/sin.vec" using 7:11 title "lx".i

set output "sin_vec.png"
replot

set title  "Nek5000 - BP1 - Scalar" font ",16"
plot for [i=start:end] "lx".i."/sin.sca" using 7:11 title "lx".i

set output "sin_scalar.png"
replot

#!/bin/sh

rm converge.pdf

gnuplot <<EOF
set terminal pdf
set output "converge.pdf"

set xlabel "Iteration Number"
set ylabel "Difference from Previous"

set xrange [0:2000]
set yrange [0:0.01]

$(cat /dev/tty)

plot $( for i in $*; do echo -n "'$i' title '$(head -n1 $i | sed -e 's/#//g')', "; done)

EOF


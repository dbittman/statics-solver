gcc main.c -o main -lm -O3
./main 600
rm out.pdf 2>/dev/null

gnuplot <<EOF
set term pdf
set output "out.pdf"
set pm3d map
splot "out" using 2:1:3:3 with image

reset

set contour
set cntrparam cubicspline
set isosamples 30,30
set pm3d implicit
splot 'out' with lines

EOF
exit 1
gnuplot <<EOF
set pm3d map
do for [c=0:1000] {
do for [i=1:10] { splot sprintf('out%d', i) using 2:1:3:3 with image; pause 0.25 }
}
EOF


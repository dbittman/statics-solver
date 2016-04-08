#!/bin/sh
VALS=()
for i in $(seq 300); do
	echo $i / 300
	VALS+=($(./solve $* | grep completed | sed -re 's/[^\(]*\(([0-9\.]+).*$/\1/p' | head -n1))
done
#echo ${VALS[*]}

echo ${VALS[*]} | awk '{ s+=$1 ; sq+=($1 * $1)} END {print "ave:",s/NR,"sq:",sq/NR,"stddev:",sqrt(sq/NR - (s/NR)*(s/NR))}' RS=" "


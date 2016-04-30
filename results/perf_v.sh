#!/bin/sh
#set -x
VALS=()
ERRVALS=()
for i in $(seq 300); do
	echo $i / 300
	RES=`./solve $*`
	SPEED=$(echo "$RES" | grep completed | sed -re 's/[^\(]*\(([0-9\.]+).*$/\1/p' | head -n1)
	ERR=$(echo "$RES" | grep 'err per cell' | sed -re 's/err per cell = ([0-9\.]+)$/\1/p' | head -n1)
	echo "$SPEED"
	echo "$ERR"
	VALS+=($SPEED)
	ERRVALS+=($ERR)
done
#echo ${VALS[*]}

echo ${VALS[*]} | awk '{ s+=$1 ; sq+=($1 * $1)} END {print "ave:",s/NR,"sq:",sq/NR,"stddev:",sqrt(sq/NR - (s/NR)*(s/NR))}' RS=" "
echo ${ERRVALS[*]} | awk '{ s+=$1 ; sq+=($1 * $1)} END {print "ave:",s/NR,"sq:",sq/NR,"stddev:",sqrt(sq/NR - (s/NR)*(s/NR))}' RS=" "


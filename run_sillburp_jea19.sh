#!/bin/bash

# USES SILLBURP TO CALCULATE GREENHOUSE GAS EMISSIONS
# FOR SILLS OF THICKNESS 50, 100, 200, 300, 400 M
# INTRUDING AT DEPTHS 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 KM

# THESE WERE THE RUNS USED IN JONES ET AL., NATURE COMMUNICATIONS, 2019


SILLBURP=/Users/smj/Dropbox/smj/work/papers/sill_papers/burp/burp_new/sillburp
CALCULATE=44

REACTIONS=7/21/55
GEOTHERM=10/0.03/1.0e-5
MAGHOST="-A1250 -D2110e-9 -W0.01/0.01"

#
# 0.5 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 11 -o $CALCULATE -eq 16 ]; then
	DEPTH=500
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<5000)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 12 -o $CALCULATE -eq 16 ]; then
	DEPTH=500
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X50/500 -J7000 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M50000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 13 -o $CALCULATE -eq 16 ]; then
	DEPTH=500
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X100/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 14 -o $CALCULATE -eq 16 ]; then
	DEPTH=500
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -J7000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/700 -O0.1 -J7000 -M250000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<250000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 15 -o $CALCULATE -eq 16 ]; then
	DEPTH=500
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/900 -O0.1 -M25000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<250000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi


#
# 1.0 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 11 -o $CALCULATE -eq 16 ]; then
	DEPTH=1000
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<5000)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 22 -o $CALCULATE -eq 26 ]; then
	DEPTH=1000
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X50/500 -J7000 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M50000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<50000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 23 -o $CALCULATE -eq 26 ]; then
	DEPTH=1000
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X100/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 24 -o $CALCULATE -eq 26 ]; then
	DEPTH=1000
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -J7000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/700 -O0.1 -J7000 -M250000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<250000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 25 -o $CALCULATE -eq 26 ]; then
	DEPTH=1000
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M100 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X20/900 -O0.1 -M10000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100 && $1<10000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>10000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi


#
# 1.5 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 31 -o $CALCULATE -eq 36 ]; then
	DEPTH=1500
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<5000)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 32 -o $CALCULATE -eq 36 ]; then
	DEPTH=1500
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M1000 > z${THICK}m.g10m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X2/300 -J7000 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 33 -o $CALCULATE -eq 36 ]; then
	DEPTH=1500
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X100/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 34 -o $CALCULATE -eq 36 ]; then
	DEPTH=1500
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -J7000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/700 -O0.1 -J7000 -M250000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 35 -o $CALCULATE -eq 36 ]; then
	DEPTH=1500
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M100 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X20/900 -O0.1 -M10000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100 && $1<10000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>10000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi


#
# 2.0 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 41 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -J2500 -O0.1 -M500 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -J2500 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<5000)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 42 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=100
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X50/500 -J7000 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M50000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 43 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X100/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 250 M THICK SILL
if [ $CALCULATE -eq 47 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=250
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X125/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X12/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 44 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/800 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X5/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 45 -o $CALCULATE -eq 46 ]; then
	DEPTH=2000
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M500 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/450 -O0.1 -M100000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>500)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi


#
# 2.5 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 51 -o $CALCULATE -eq 56 ]; then
	DEPTH=2500
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 52 -o $CALCULATE -eq 56 ]; then
	DEPTH=2500
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M1000 > z${THICK}m.g10m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X2/300 -J7000 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 53 -o $CALCULATE -eq 56 ]; then
	DEPTH=2500
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X100/1500 -O0.1 -J8000 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/1000 -O0.1 -J8000 -M100000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 54 -o $CALCULATE -eq 56 ]; then
	DEPTH=2500
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M500 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/450 -O0.1 -M100000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>500 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100000 && $1<1000000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 55 -o $CALCULATE -eq 56 ]; then
	DEPTH=2500
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M500 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X10/450 -O0.1 -M100000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>500 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100000 && $1<1000000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi


#
# 3.0 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 61 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 62 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X50/500 -J7000 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M50000 > z${THICK}m.g10m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<50000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 63 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X125/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X12/700 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X4/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 250 M THICK SILL
if [ $CALCULATE -eq 67 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=250
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X125/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X12/700 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X4/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 64 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/800 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X5/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 65 -o $CALCULATE -eq 66 ]; then
	DEPTH=3000
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X20/900 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100 && $1<50000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi



#
# 3.5 KM DEEP INTRUSION 
#

# 50 M THICK SILL
if [ $CALCULATE -eq 71 -o $CALCULATE -eq 76 ]; then
	DEPTH=3500
	THICK=50
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X25/750 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/900 -O0.1 -M5000 > z${THICK}m.g5m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g5m.time >> z${THICK}m.dec.time
	date
fi

# 100 M THICK SILL
if [ $CALCULATE -eq 72 -o $CALCULATE -eq 76 ]; then
	DEPTH=3500
	THICK=100
	date
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X5/700 -J7000 -O0.1 -M1000 > z${THICK}m.g10m.time
	$SILLBURP -Z100/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -N$REACTIONS -X2/300 -J7000 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250 && $1<100000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	date
fi

# 200 M THICK SILL
if [ $CALCULATE -eq 73 -o $CALCULATE -eq 76 ]; then
	DEPTH=3500
	THICK=200
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X125/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X12/700 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X4/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 300 M THICK SILL
if [ $CALCULATE -eq 74 -o $CALCULATE -eq 76 ]; then
	DEPTH=3500
	THICK=300
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X150/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X15/800 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X5/140 -O0.1 -M1000000 > z${THICK}m.g30m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>250)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g30m.time >> z${THICK}m.dec.time
	date
fi

# 400 M THICK SILL
if [ $CALCULATE -eq 75 -o $CALCULATE -eq 76 ]; then
	DEPTH=3500
	THICK=400
	date
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X200/1500 -O0.1 -M250 > z${THICK}m.g1m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X20/900 -O0.1 -M50000 > z${THICK}m.g10m.time
	$SILLBURP -Z${THICK}/$DEPTH $MAGHOST -N$REACTIONS -g$GEOTHERM -X2/600 -O0.1 -M1000000 > z${THICK}m.g50m.time
	cat z${THICK}m.g1m.time > z${THICK}m.dec.time
	awk '{if(NR>1 && $1>100 && $1<50000)print $0}' z${THICK}m.g10m.time >> z${THICK}m.dec.time
	awk '{if(NR>1 && $1>50000)print $0}' z${THICK}m.g50m.time >> z${THICK}m.dec.time
	date
fi



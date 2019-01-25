#!/bin/bash

# input string
inputstr=$1
odir=$2

# parse input string
file=${inputstr##*/}
baseid=${file%%.dat}
title=${baseid%*_N*}
cfgseed=${baseid#*seed*}
str1=${baseid%*_seed*}
N=${str1#*N*}

echo file = $file
echo baseid = $baseid
echo title = $title
echo str1 = $str1
echo N = $N
echo cfgseed = $cfgseed

# prepare xyz and cfg string
xyzstr=$odir"$title"_N"$N"_cfg"$cfgseed"_vseed
cfgstr=$xyzstr


# run new vseed, print xyz and cfg
NT=200000
T0=0.0001
vseed=1

# finish output strings
xyzstr=$xyzstr$vseed.xyz
cfgstr=$cfgstr$vseed.dat

# code directories
srcdir=~/_pv/sim/generate_packings/src
maindir=~/_pv/sim/generate_packings/main
binf=res_md.o

# compile 
g++ -I $srcdir $maindir/res_rcp_trajectory.cpp $srcdir/*.cpp -o $binf

if [[ ! -f $binf ]]
then
    echo bin file does not exist, must have been a compilation error....
    echo ending...
    exit 1
fi

# run
echo running binary $binf
./$binf $N $NT $T0 $vseed $inputstr $xyzstr $cfgstr

# remove binary
rm -f $binf

# ending
echo finished running, hopefully no errors







#!/bin/bash

# ==================================
# script to run residue_vel_run.cpp
# ==================================
echo Running code to pack a systems of residues, run NVE MD and output velocities

# IMPORTANT: MUST DEFINE DIRECTORY WITH generate_packings GIT FOLDERS
workingdir=~/_pv/sim/generate_packings

# code directories (DO NOT EDIT unless src code or main file is movied)
srcdir="$workingdir"/src
maindir="$workingdir"/examples/res_vacf
mainf=residue_vel_run.cpp

# EDIT THIS TO DEFINE OUTPUT DIRECTORY
outdir="$workingdir"

# INPUT PARAMETERS
N=$1					# Number of residues
phi0=$2					# initial packing fraction (for compression)
dphi=$3					# packing fraction step (for compression)
dphiDM=$4				# packing fraction increase (for MD run after compression)
vacfNT=$5				# number of MD steps (after compression)
vacfT0=$6				# temperature for MD run (after compression)
tsave=$7				# number of steps to skip between saves in MD run (after compression)
input_str=$8			# full path to input file

# TURN NT and tsave into powers of 2
let vacfNT=2**"$vacfNT"
let tsave=2**"$tsave"

# determine input seed
file_name=${input_str##*/}
file_base=${file_name%%.dat}
seed=${file_base#*seed*}


# CREATE OUTPUT FILE NAMES
echo -- creating output files

# jammed configuration file
cfg_str="$outdir"/res_rcp_config_N"$N"_seed"$seed".dat

# jammed configuration state file
stat_str="$outdir"/res_rcp_stat_N"$N"_seed"$seed".dat

# MD trajectory file
vel_str="$outdir"/res_rcp_vel_N"$N"_NT"$vacfNT"_T0"$vacfT0"_tsave"$tsave"_dphi"$dphiDM"_seed"$seed".dat

# echo file strings
echo -- config file" 		: $cfg_str"
echo -- stat file" 			: $stat_str"
echo -- velocity file"		: $vel_str"

# BINARY FILE
binf="res_vel.o"
rm -f "$binf"

# COMPILE CODE
echo -- compiling...
g++ --std=c++11 -I "$srcdir" "$maindir"/"$mainf" "$srcdir"/*.cpp -o "$binf"


# CHECK COMPILATION
if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compile failed.
    exit 1
else
	echo -- cmpilation successful!
fi

# run code
echo Running!
echo ..
echo ..
param_str="$N $phi0 $dphi $dphiDM $vacfNT $vacfT0 $tsave $seed $input_str $cfg_str $stat_str $vel_str"
echo ./$binf $param_str
echo ..
echo ..
./$binf $param_str
echo code completed!





# ====================
#       INPUTS
# ====================
# 1. N
# 2. phi0
# 3. dphi
# 4. dphiDM
# 5. vacfNT
# 6. vacfTO
# 7. tsave
# 8. input string





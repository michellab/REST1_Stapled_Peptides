#! /bin/bash

#SBATCH --job-name=checksoluteT
#SBATCH -p GTX
#SBATCH --gres=gpu:1
#SBATCH --time 48:00:00

module load fftw/3.3.4
module load amber/16
module load cuda/7.5

source  /home/common/AMBER16/amber16/amber.sh 
export LD_LIBRARY_PATH=/users/common/miniconda2/lib/:$LD_LIBRARY_PATH





if [ "$1" != "" ]; then
    inputfile=$1
    echo 'running file' $1 
else
    echo "Positional parameter 1 is empty"
fi

if [ "$2" != "" ]; then
    temp=$2
    echo 'running at ' $2 
else
    echo "Positional parameter 2 is empty (temperature)"
fi


cp ~/submitambersim/equilibrateT.in .
cp ~/submitambersim/equilibrate.in equilibrate1.in 
cp ~/submitambersim/equilibrate2.in  .
cp ~/submitambersim/pre-equilibrate2.in  .

cp ~/submitambersim/prodT-noexg.in prodT.in

sed -i s/\$TEMP/$temp/g   prodT.in
sed -i s/\$TEMP/$temp/  equilibrateT.in


pmemd.cuda -O -i ~/submitambersim/min-restrain.in -p   $inputfile\.prmtop  -ref  $inputfile\.rst7 -c $inputfile\.rst7 -r min.nc;
pmemd.cuda -O -i  equilibrate1.in -p $inputfile\.prmtop  -c  min.nc  -r eq.nc  # 0 -->> 300 K   
sander -O -i pre-equilibrate2.in -p $inputfile\.prmtop  -c eq.nc -r pre-eq2.nc     #  NPT 
pmemd.cuda  -O -i equilibrate2.in -p $inputfile\.prmtop  -c pre-eq2.nc -r eq2.nc     #  NPT
pmemd.cuda -O -i  equilibrateT.in -p $inputfile\.prmtop  -c  eq2.nc  -r eq3.nc   ## 300 -->> target_temp
pmemd.cuda  -O -i prodT.in -p $inputfile\.prmtop  -c eq3.nc -r restart.nc -x traj.nc
python ~/submitambersim/prepare_AMD.py $inputfile\.prmtop $inputfile\.rst7
sbatch ~/submitambersim/runAMBERs9aMD.bsh   $inputfile  $temp

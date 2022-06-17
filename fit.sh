#!/bin/sh

if(($#<5 || $#>9)); then
    echo "ERROR: wrong input."
    echo "Usage: bash ${0} <mesh-factor> <n.atoms> <qmin> <qnmax> <dq> [gamma] [t1] [t2] [fitting-mode]"
    echo ''
    echo "This script calculates fits Re[F(q,t)] using the script fit.py for the q's in the given interval [qnmin,qnmax] every dq; qmin,qmax are positive integers in units of the mesh grid pi/L*meshfactor. The remaining arguments are the same as for fit.py."
    echo ''
    exit 1
fi

M=$1
NAT=$2
qmin=$3
qmax=$4
dq=$5
args=""$6" "$7" "$8" "$9 #pass extra arguments (may be empty) separated with spaced

for ((q=$qmin;q<=$qmax;q+=dq)); do
    printf -v phifile "phi_q%03d_m%03d_nat%d.dat" $q $M $NAT
    echo ""
    printf "q%03d ...\n" $q
    echo ""
    python3 fit.py $phifile $args
done

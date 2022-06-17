#!/bin/sh

if(($#<3 || $#>5)); then
    echo "ERROR: wrong input."
    echo "Usage: bash ${0} <mesh-factor> <n.atoms> <qmin> [qnmax=qnmin] [dq=1]"
    echo ''
    echo "This script calculates the Re[F(q,t)] and Im[F(q,t)] for the q's in the given interval [qnmin,qnmax]; qmin,qmax are positive integers in units of the mesh grid pi/L*meshfactor (meshfactor is a positive integer). n.atoms is the grouped atoms for the center of mass of the molecules. Average is over q-vector direction and initial time t0."
    echo ''
    exit 1
fi

M=$1
NAT=$2
qmin=$3
qmax=qmin
dq=1
if(($#>=4)); then
    qmax=$4
    if(($#==5)); then
        dq=$5
    fi
fi

gcc -Wall -O3 average.c -lm -o average.x

for ((i=$qmin;i<=$qmax;i+=dq)); do
    printf -v logdata "fqt_q%03d_m%03d_nat%d-log.dat" $i $M $NAT
    printf -v lindata "fqt_q%03d_m%03d_nat%d-lin.dat" $i $M $NAT
    printf -v infile "fqt_q%03d_m%03d_nat%d.dat" $i $M $NAT
    printf -v outfile "fqt_ave_q%03d_m%03d_nat%d.dat" $i $M $NAT
    printf -v phifile "phi_q%03d_m%03d_nat%d.dat" $i $M $NAT
    printf "q%03d ...\n" $i
    
    cat $logdata $lindata > $infile
    #if(($?==1)); then break; fi
    # sort the data by Dt and send it to awk script for the average
    #LC_ALL=C sort -g -k 2 $infile | gawk -f average.awk > $outfile
    LC_ALL=C sort -g -t " " -k 2 $infile > temp
    #echo "n. of rows with dt-idx==0:"
    #gawk '{if($2==0) x++;}END{print x}' temp
    Nlines=$(cat temp | wc -l)
    ./average.x temp $Nlines > $outfile
    gawk '{print $1,$7,$8}' $outfile > $phifile
    rm $infile
#    python3 plot_fqt.py $outfile
done

rm temp
rm average.x

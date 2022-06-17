# list configurations | order by time | remove .tpr, .gz, ... | you get a list of configuration filenames
#ls ../CnfLT500*-* | sort -n -t "-" -k 2 | awk -F "." '(NR>1){if($4=="") print $0}' > tmp
rm s*-vs-t*dat

M=1
qmax=30

for NAT in {1,3}
do
	printf -v pattern "struct-????_q%03d_m%03d_nat%d.dat" $qmax $M $NAT
	ls $pattern | sort -n -t "q" -k 1 | while read nome
	do
	    #get time index without zeros in front
	    timeidx=$(echo $nome | awk -F "q" '{print $1}' | awk -F "-" '{x=int($2); print x}')
	    #get time (ps) and box size (nm)
	    nrow=$((1+$timeidx))
	    idx=$(head -n $nrow listatempi.dat | tail -n 1 | awk '{print $1}')
	    t=$(head -n $nrow listatempi.dat | tail -n 1 | awk '{print $2}')
	#    conf=../Cnf*-${idx}
	#    t=$(head -n 1 $conf | awk -F "t=" '{print $2}' | awk -F "step=" '{x=$1; print x}')
	#    L=$(tail -n 1 $conf | awk '{x=$1; print x}')
	#    echo "{print \$1*3.141592/$L * $M}" > awkfile
	    cat $nome | head -n 8 | while read qn sq sqerr
	    do
	#        q=$( echo $qn | awk -f awkfile )
	#	echo $q
		echo $t $sq $sqerr >> s${qn}-vs-t_nat${NAT}.dat
	    done
	done
done

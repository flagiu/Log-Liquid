if(($#!=1))
then
	echo "Usage $0 <mesh-factor>"
	echo ""
	exit 1
fi

M=$1
echo "s/M = [0-9]^*/M = "$M"/" > script
sed -i -f script avestruct.awk #replace mesh factor in avestruct.awk
rm script
for NAT in {1,3}
do
	for qn in {30,462}
	do
		printf -v pattern "struct-????_q%03d_m%03d_nat%d.dat" $qn $M $NAT
		printf -v outfile "sq_q%03d_m%03d_nat%d.dat" $qn $M $NAT
		ls $pattern >> /dev/null
		if(($?==0)) #check if there are files for that pattern
		then
			awk -f avestruct.awk $pattern > $outfile
			echo ""$outfile" created"
		fi
	done
	# now collect all different qn's
	printf -v pattern "sq_q???_m%03d_nat%d.dat" $M $NAT
	printf -v outfile "sq_m%03d_nat%d.dat" $M $NAT
	ls $pattern >> /dev/null
	if(($?==0)) #check if there are files for that pattern
	then
		cat $pattern > $outfile
		echo ""$outfile" created"
	fi
done

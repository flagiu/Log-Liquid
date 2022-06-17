for varname in {"fq","tauA","beta","tau","omega"}; do rm ${varname}.dat; done
ls fit_q*.log | sort -n | while read nome
do
    q=$( echo $nome | gawk -F "_q" '{print $2}' | gawk -F ".log" '{print $1}' )
    for varname in {"fq","tauA","beta","tau","omega"}; do
	var=$( grep ${varname}"=" $nome | gawk -F "=" '{print $2}' | gawk -F "+-" '{print $1,$2}' )
	echo $q $var >> ${varname}.dat
    done
done

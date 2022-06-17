echo "# t(ps) : V(nm^3)" > V.dat
cat lista.dat | while read file idx t
do
	V=$(tail -n 1 $file | awk '{print $1*$2*$3}')
	echo $t $V >> V.dat
done
NPC=$( head -n 1 ../generatempi.f | awk -F 'npc=' '{print $2}' | awk -F ',' '{print $1}')
echo "( (NR % (${NPC}/8)) ==0 ){N++;V+=\$2;sV+=\$2*\$2}END{print V/N, sqrt( (sV/N - (V/N)**2)/N )}" > awkfile
echo "# V avg (nm^3) : V avg error" > Vavg
awk -f awkfile V.dat >> Vavg
rm awkfile

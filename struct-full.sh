if(($#!=1)); then
	echo "Usage: $0 <t0-idx, starting from 0>"
	exit 1
fi
# il primo picco è a q=200
gcc -O3 -D MESHFACT=1 struct.c -lm -o struct.x
./struct.x 10 180 10 $1
./struct.x 181 200 1 $1
./struct.x 201 219 1 $1
./struct.x 220 300 10 $1
# il secondo picco è a q=350=2*175
gcc -O3 -D MESHFACT=2 struct.c -lm -o struct.x
./struct.x 155 200 1 $1
./struct.x 205 300 5 $1

printf -v root "struct-%04d" $1
cat ${root}q180.dat  ${root}q200.dat  ${root}q219.dat  ${root}q300.dat ${root}q400.dat ${root}q600.dat >  ${root}.dat
rm struct-*q*.dat
echo "Structure factor for t0=$1 completed"
awk -f avestruct.awk struct-????.dat > sq.dat
#xmgrace -settype xydy struct-$1.dat &

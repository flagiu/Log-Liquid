if(($#!=1))
then
    echo "Usage: $0 <qvector file>"
    echo ""
    echo "The program outputs the min. and max. modulus wavenumber found in the given file."
    echo ""
    exit 1
fi

file=$1
# First loop: find the lower wavenumber
res=1
for((n=0;res!=0;n++))
do
    grep $n $file >> /dev/null
    res=$?
    #echo "n="$n"  res="$res #debugging
done
#the error disappeared 1 step ago
n=$((n-1))
echo "Min modulus wavenumber is "$n
#Second loop: find the highest wn
for((;res==0;n++))
do
    grep $n $file >> /dev/null
    res=$?
    #echo "n="$n"  res="$res #debugging
done
#the error occurred 1 step ago
#so the last non-error is 2 steps earlier 
n=$((n-2))
echo "Max modulus wavenumber is "$n

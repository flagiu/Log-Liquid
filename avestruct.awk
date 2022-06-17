# media e dev.standard file di 2 colonne del tipo x vs.y
# i file possono avere n.righe diverse, ma devono essere allineate per x
# es: il primo ha x=1,2,...,10; il secondo estende x=1,2,...,13
BEGIN {
    #ncols=2
    pi = atan2(0, -1)
    L = 43.32963 #nm # PROBLEM FOR NPT ENSEMBLE... L=L(t)
    M = 1   #mesh factor (replaced automatically in avestruct.awk)
    mesh=pi/L * M
}
{
    nfiles[FNR] ++
    sum[FNR,1] = $1*mesh
    sum[FNR,2] += $2
    sum[FNR,3] += $2*$2
    #for (i = 2; i <= ncols; i++) sum[FNR,i] += $i
    if (FNR > maxnr) maxnr = FNR
}
END {
    for (line = 1; line <= maxnr; line++)
    {
	printf " %f", sum[line,1];
	#for (col = 2; col <= ncols; col++)
	#    printf "  %f", sum[line,col]/nfiles;
	ave = sum[line,2]/nfiles[line];
	std = sqrt(sum[line,3]/nfiles[line] - ave*ave);
	if(nfiles[line]>1) std/=sqrt(nfiles[line]-1);
	printf " %f %f\n", ave, std;
    }
}

# Writes a list of filename,step,time into lista.dat
print(" ")
with open("listaconf.dat","r") as f1:
    lines = f1.readlines()
    if len(lines)==0:
        print("listaconf.dat is empty")
    
    fout=open("lista.dat","w")
    ftempi=open("listatempi.dat","w")
    for line in lines:
        fname=line.strip('\n')
        print("\rprocessing %s ... "%fname,end='')
        try:
            step = int(''.join(filter(str.isdigit,fname.split('-')[-1])))
            with open(fname,"r") as f2:
                firstline = f2.readline()
                t=float(firstline.split('t=')[-1].split('step=')[0])
                fout.write("%s %d %f\n"%(fname,step,t))
                ftempi.write("%d %f\n"%(step,t))
        except ValueError:
            print(" Skipping invalid name")
    fout.close()
    ftempi.close()
print("done\n")

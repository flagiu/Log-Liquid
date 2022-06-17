#include "lib.h"

//variable qnmax is an integer defined in lib.h

void printError(char *myname);
void restart(char*, int*, int*, int*);
int getTimes(int*, double*);
int getQvectors(char*, int[MAXQARR][3]);
void savelog(char *, int, int, int);
void updateTerminal(double,int,int, double,int,int, int[3],int,int);
void get_exp(double *exp, int step, int qnmax);

#define ALL 0 //t0 mode
#define SINGLE 1

void printError(char *myname) {
  fprintf(stderr,
	  "\nNormal Usage: %s <qmin> <qmax> <dq> \n\nQuesto programma calcola lo Static Structure Factor S(q)=F(q,Dt=0) a partire dai frame corrispondenti al pattern ../Cnf*-*. Per tutti i vettori di modulo q (vedi cartella QVECTORS) La F((qx,qy,qz),t0,Dt=0) viene appesa al file struct-XXXX_qXXX.dat. Questi dati dovranno essere mediati esternamente su t0.\n  q-number is an integer in mesh units (>=2)\n"
	  "\nSingle-t0 mode: %s <qmin> <qmax> <dq> <t0-idx>\n\n Se vengono dati 4 argomenti anziché 3, il programma calcola solo S(q) al tempo t0, corrispondente alla configurazione distante t0-idx configurazioni da quella iniziale.\n"
	  "\nNota: se qmax<=%d il programma precalcola tutti gli esponenziali, altrimenti li calcola sul momento.\n\n", myname,myname,MAXNMAX);
  exit(1);
}

int main(int argc, char **argv) {
  int j,l,jj,i,Nq, nmax,nmin,dn,n, steps[MAXCNF], Ncnf, ncycles, qarr[MAXQARR][3];
  int t0restart=0,Qvecrestart=0, lenT0cycle,stepT0cycle;
  int first=1, Nqvec, printexp=0, t0mode, onfly;
  double ts[MAXCNF], *e0, *e1, rho0[2], rhot[2], fqt[2], *S, *S2, qn;
  char logfile[15], outfile[50], qvecfile[100];
  FILE *fp;
  if (argc!=4 && argc!=5) printError(argv[0]);
  Ncnf = getTimes(steps, ts); //returns n. of configurations and their times (ps), and initialized cnfroot, N, L !!!!
  ncycles=Ncnf/NPC;
  nmin=atoi(argv[1]);
  nmax=atoi(argv[2]);
  dn=atoi(argv[3]);
  sprintf(logfile, "struct.logc");
  if(argc==4)
    {
      t0mode=ALL;
      fprintf(stderr,"t0-mode: ALL\n");
      restart(logfile, &t0restart,NULL, &Qvecrestart);
    }
  else
    {
      t0mode=SINGLE;
      fprintf(stderr,"t0-mode: SINGLE\n");
      t0restart=atoi(argv[4]);
    }
  if (nmax<=0 || nmin<=0 || dn<=0 || nmin>nmax || t0restart<0) printError(argv[0]);
  lenT0cycle  = t0mode==ALL ? ncycles : 1;
  stepT0cycle = t0mode==ALL ?     NPC : 1;
  Nq=(nmax-nmin)/dn+1;
  onfly=(nmax>MAXNMAX);
  qnmax=floor(DN*nmax);
  if(onfly)
    e0=(double *)malloc(3*NMOL*sizeof(double)); //store only r
  else
    e0=(double *)malloc(EXPN*sizeof(double)); //store exp(Dq*n*r)
  S  = (double *)calloc(Nq,sizeof(double));
  S2 = (double *)calloc(Nq,sizeof(double));
  fprintf(stderr, "\nStarting..\n");
  for (jj=0; jj<lenT0cycle; jj++) {
    j=t0restart+jj*stepT0cycle;
    if(onfly) get_coords(e0, steps[j]);
    else get_exp(e0, steps[j], qnmax);
    fprintf(stderr,"t0=%f (%d/%d)\n",ts[j],jj+1,lenT0cycle);
    for(i=0;i<Nq;i++){
      n=nmin+dn*i;
      //qn=floor(DN*n);
      sprintf(qvecfile, "QVECTORS/qvector.%03d",n);
      Nqvec = getQvectors(qvecfile, qarr);
      for (l=0; l<Nqvec; l++) {
	calcRHORHO(fqt, e0,e0, qarr[l], rho0,rhot, printexp, onfly);
	S[i]+=fqt[0];
	S2[i]+=(fqt[0]*fqt[0]);
	if (first) first=t0restart=Qvecrestart=0;
      }
      S[i]/=(double)Nqvec; //mean
      S2[i]/=(double)Nqvec;
      S2[i]=sqrt((S2[i]-S[i]*S[i])/(double)Nqvec); //std.dev. of mean
      
      fprintf(stderr," qn=%4d (%4d/%4d); S(q) = %.5g +/- %.5g\n",n*MESHFACT,i+1,Nq,S[i],S2[i]);
    }
    sprintf(outfile, "struct-%04d_q%03d_m%03d_nat%d.dat",j,nmax,MESHFACT,NAT);
    fp=fopen(outfile, "w");
    for(i=0;i<Nq;i++)
      {
	fprintf(fp,"%d %f %f\n",(nmin+dn*i)*MESHFACT,S[i],S2[i]);
	S[i]=S2[i]=0.0; //reset averages
      }
    fprintf(stderr,"Data saved on %s\n",outfile);
    fclose(fp);
  }
  
  free(e0);
  free(e1);
  free(S);
  free(S2);
  return 0;
}

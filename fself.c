#include "lib.h"

//variable qnmax is an integer defined in lib.h

void printErr(char *myname);
void restart(char*, int*, int*, int*);
int getTimes(int*, double*);
int getQvectors(char*, int[MAXQARR][3]);
void savelog(char *, int, int, int);
void updateTerminal(double,int,int, double,int,int, int[3],int,int);
void get_exp(double *exp, int step, int qnmax);

int main(int argc, char **argv) {
  int j,i,l,jj,ii, nmax,nmax1,nmax2,dnmax,n, steps[MAXCNF], Dtmode, Ncnf, ncycles, qarr[MAXQARR][3];
  int t0restart=0,Dtrestart=0, *restarts, Nnmax, *checkptReached;
  int lenT0cycle,stepT0cycle, lenDTcycle,stepDTcycle;
  int first=1, Nqvec, printexp=0, onfly;
  double ts[MAXCNF], Dt, *e0, *e1, Fself[2];
  char logfile[50], Dtmodestr[4], outfile[50], qvecfile[100];
  FILE *fp;
  if (argc<3 || argc>5) printErr(argv[0]);
  Ncnf = getTimes(steps, ts); //n. of configurations available
  //////// debugging mode /////////////
  /*
  if (argc == 5) {
    printexp=1;
    nx=atoi(argv[1]);
    ny=atoi(argv[2]);
    nz=atoi(argv[3]);
    qarr[0][0]=nx;
    qarr[0][1]=ny;
    qarr[0][2]=nz;
    i=atoi(argv[4]);
    qnmax= (abs(nx)>abs(ny) ? abs(nx) : abs(ny));
    qnmax= (qnmax>abs(nz) ? qnmax : abs(nz));
    e0 = (double *)malloc(EXPN*sizeof(double)); 
    e1 = (double *)malloc(EXPN*sizeof(double));
    get_exp(e0,steps[0],qnmax);
    get_exp(e1,steps[0+i],qnmax);
    calcRHORHO(fqt, e0,e1, qarr[0], rho0,rhot, printexp, 0);
    return 0;
    }*/
  ////// normal mode ////////
  Dtmode=atoi(argv[1]);
  nmax1=atoi(argv[2]);
  if (nmax1<=0) printErr(argv[0]);
  nmax2=nmax1;
  dnmax=1;
  if(argc>3) {nmax2=atoi(argv[3]); if(nmax2<nmax1) printErr(argv[0]);}
  if(argc>4) {dnmax=atoi(argv[4]); if(dnmax<1) printErr(argv[0]);}
  qnmax=floor(DN*nmax2);
  onfly=(2*nmax2>MAXNMAX);
  if (Dtmode!=LOG && Dtmode!=LIN) printErr(argv[0]);
  strcpy(Dtmodestr, (Dtmode==LOG)?"log":"lin");

  if(onfly)
    {
      e0 = (double *)malloc(3*NMOL*sizeof(double)); 
      e1 = (double *)malloc(3*NMOL*sizeof(double));
    }
  else
    {
      e0 = (double *)malloc(EXPN*sizeof(double)); 
      e1 = (double *)malloc(EXPN*sizeof(double));
    }
  ncycles= Ncnf/NPC; //n. of log time cycles
  Nnmax=(nmax2-nmax1)/dnmax+1;
  restarts = (int *)malloc(3*Nnmax*sizeof(int));
  checkptReached = (int *)calloc(2*Nnmax, sizeof(int));

  fprintf(stderr,"\nLoading checkpoints...\n");
  for(i=0;i<Nnmax;i++) {
    nmax=nmax1+i*dnmax;
    sprintf(logfile,       "q%03d_m%03d_nat%d-%s_fself.logc", nmax,MESHFACT,NAT,Dtmodestr);
    sprintf(outfile, "fself_q%03d_m%03d_nat%d-%s.dat",        nmax,MESHFACT,NAT,Dtmodestr);
    sprintf(qvecfile, "QVECTORS/qvector.%03d",nmax);
    restart(logfile, restarts+3*i+0, restarts+3*i+1, restarts+3*i+2);
    Nqvec = getQvectors(qvecfile, NULL); //just count lines
    if(i==0){
      t0restart=restarts[3*i+0];
      Dtrestart=restarts[3*i+1];
    }
    if(restarts[3*i+0]<=t0restart) {
      t0restart=restarts[3*i+0];
      if(restarts[3*i+1]<=Dtrestart) {
	Dtrestart=restarts[3*i+1];
	n=i;
      }
    }
  }
  fprintf(stderr,"...done.\nEarliest starting point is t0=%d, Dt=%d for q=%d\n\n",t0restart,Dtrestart,(nmax1+n*dnmax)*MESHFACT);
  
  lenT0cycle  = Dtmode==LOG ? ncycles : ncycles;//Ncnf-NPC;
  stepT0cycle = Dtmode==LOG ?     NPC : NPC;//1;
  for (jj=t0restart; jj<lenT0cycle; jj++) {
    fprintf(stderr,"\r[%d%% t0] ",(int)((100.0*jj)/lenT0cycle));
    j=jj*stepT0cycle;
    if(onfly) get_coords(e0, steps[j]);
    else get_exp(e0, steps[j], qnmax);
    fprintf(stderr,"\n");
    
    lenDTcycle  = Dtmode==LOG ? NPC : ncycles-1-j/NPC;
    stepDTcycle = Dtmode==LOG ?   1 : NPC;
    if(Dtmode==LIN && Dtrestart==0) Dtrestart=1; //include Dt=0 only in LOG mode
    for (ii=Dtrestart; ii<lenDTcycle; ii++) {
      i=ii*stepDTcycle;
      Dt=ts[j+i]-ts[j];
      fprintf(stderr,"\r");
      if(onfly) get_coords(e1, steps[j+i]);
      else get_exp(e1, steps[j+i], qnmax);
      fprintf(stderr," Computing Fself at t0=%f (%d/%d); Dt=%f (%d/%d) for all q's... ",ts[j],jj+1,lenT0cycle,Dt,ii+1,lenDTcycle);
      for(n=0;n<Nnmax;n++) {
	if(jj<restarts[3*n+0] || ii<restarts[3*n+1])
	  //vero finché non raggiunge t0,Dt corrispondente al checkpoint di q...
	  continue;
	if(!checkptReached[2*n+0]) {
	  checkptReached[2*n+0]=1;
	  restarts[3*n+0]=restarts[3*n+1]=0; //...poi sarà sempre falso
	}
	nmax=nmax1+n*dnmax;
	sprintf(logfile,       "q%03d_m%03d_nat%d-%s_fself.logc", nmax,MESHFACT,NAT,Dtmodestr);
	sprintf(outfile, "fself_q%03d_m%03d_nat%d-%s.dat",        nmax,MESHFACT,NAT,Dtmodestr);
	sprintf(qvecfile, "QVECTORS/qvector.%03d",nmax);
	Nqvec = getQvectors(qvecfile, qarr);
	fp=fopen(outfile, "a");
	for (l=0; l<Nqvec; l++) {
	  if(l<restarts[3*n+2]) continue;//vero finché non raggiunge il q-vector corrispondente al checkpoint di q...
	  if(!checkptReached[2*n+1]) {
	    checkptReached[2*n+1]=1;
	    restarts[3*n+2]=0; //...poi sarà sempre falso
	  }
	  calcFself(Fself, e0,e1, qarr[l], onfly);
	  savelog(logfile, jj,ii,l);
	  fprintf(fp, "%d %d %18e %d %d %d %18e %18e\n",
		  j,i,Dt,qarr[l][0],qarr[l][1],qarr[l][2],Fself[0],Fself[1]
		  );
	  //updateTerminal(ts[j],jj+1,lenT0cycle,Dt,ii+1,lenDTcycle,qarr[l],l+1,Nqvec); //troppo veloce, questo output è inutile
	  if (first) first=t0restart=Dtrestart=0;
	}
	fflush(fp);
	fclose(fp);
      }
    }
    fprintf(stderr,"done\n");
  }
  
  free(e0);
  free(e1);
  free(restarts);
  return 0;
}


void printErr(char *myname) {
  fprintf(stderr,
	  "\nNormal Usage: %s <Dt-mode> <q-min> [q-max=q-min] [dq=1]\n\nQuesto programma calcola la Self Intermediate Scattering Function Fs(q,Dt) a partire dai frame corrispondenti al pattern ../Cnf*-*. Per tutti i vettori di modulo q (vedi cartella QVECTORS), per q compreso in [qmin,qmin+dq,...,qmax]. La Fs((qx,qy,qz),t0,Dt) viene appesa al file fself_qXXX.dat. Questi dati dovranno essere mediati esternamente su t0 e sul vettore q.\n Dt-mode\t %d (log) per dt piccoli; %d (lin) per dt grandi.\n"
	  " -Dt-mode: 0 logarithmic times, 1 linear times;\n"
	  " -q-min,q-max,dq: integers in units of MESH*dq;\n\n"
	  "#(INACTIVE!) Debugging mode: %s <qnx> <qny> <qnz> <Dt-idx>\n"
	  "#\n"
	  "# Se vengono dati 4 argomenti anziché 2, il programma calcola solo la Fs(2pi/L*(qnx,qny,qnz),0,Dt) con Dt corrispondente alla configurazione distante Dt-idx configurazioni da quella iniziale. Stampa sullo stdout i valori di exp(iqx*rx) per ogni particella, per qx!=0 (è pensato per vettori multipli di 100,010,001).\n\n", myname,LOG,LIN,myname);
  exit(1);
}

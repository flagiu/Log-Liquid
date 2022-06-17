#include"lib.h"

void printErr(char *myname) {
  fprintf(stderr, "\nUsage: %s [t0] [t1] \n\nQuesto programma calcola il MSD(t) a partire dai frame corrispondenti contenuti nel file lista.dat, a MULTIPLI DI NPC.\n t0\ttempo iniziale (opzionale), in ps\n t1\ttempo finale (opzionale), in ps\n\n", myname);
  exit(1);
}

int main(int argc, char **argv) {
  int a,i,i0,i1,j,k,xi,nt,Ncnf,idx,idx_old, steps[MAXCNF], ncycles, *Navg, maxncycles, ti;
  float rr[3], r1[3], d, foo[3];
  double dt,t1,ts[MAXCNF], *msd,*msd2,dr,rc[3], *r, dr2,dr4, tcycle;
  char str[10],cnffile[50];
  char outfile[30];
  FILE *fcnf,*fout;
  Ncnf = getTimes(steps, ts); //return n. of configurations, their time (ps), and initializes cnfroot, N, L !!!!!!!!
  i0=0;
  i1=Ncnf-1;
  fprintf(stderr,"Available time window: [%lf, %lf] ps\n",ts[i0],ts[i1]);
  if(Ncnf>NPC) tcycle=ts[NPC-1]-ts[0];
  else {
    fprintf(stderr,"Error: less than 1 log-cycle available. Exiting\n");
    exit(1);
  }
  if(argc>3) printErr(argv[0]);
  else if (argc>=2) {
    t1=atof(argv[1]);
    while(ts[i0]<t1 && i0<Ncnf-1) i0++;
    i0=(i0/NPC)*NPC; //prendi il tempo zero piu' vicino da sotto
    if (argc==3) {
      t1=atof(argv[2]);
      while(ts[i1]>t1 && i1>0) i1--;
      i1=((i1-i0)/NPC)*NPC+i0;
    }
  }
  ncycles=(i1-i0)/NPC;
  fprintf(stderr,"You chose: \t\t[%f, %f] ps\n => %d cycles (%f ps each)\n",ts[i0],ts[i1], ncycles,tcycle);
  /* Account for non-continuous-in-time cycles (e.g. 2nd cycle was skipped) */
  maxncycles=floor((ts[i1]-ts[i0])/tcycle)+1;
  fprintf(stderr,"Maximum Dt in the selected interval:\n maxncycles=%d => maxDt=%f ps\n",maxncycles,(maxncycles-1)*tcycle);
  
  //restart(logfile, &trestart);
  r=(double *)malloc(NMOL*3*ncycles*sizeof(double)); // trajectory of each particle
  msd=(double *)calloc(maxncycles,sizeof(double)); //maxncycles+1 !!!
  msd2=(double *)calloc(maxncycles,sizeof(double));
  Navg=(int *)calloc(maxncycles,sizeof(int));
  fprintf(stderr, "Allocated maxncycles slots for MSD.\n");

  fprintf(stderr,"\nReconstrucing trajectory and MSD: \n");
  // 1) LOOP ON t1
  for(i=0;i<ncycles;i++) {
    t1=ts[i0+i*NPC];
    sprintf(str,"-%d",steps[i0+i*NPC]);
    strcpy(cnffile,cnfroot);
    strcat(cnffile,str);
    //fprintf(stderr,"\nopening file %s\n",cnffile);
    fcnf=fopen(cnffile,"r");
    read_N_L(fcnf, &N, L);
    // 2) LOOP ON MOLECULES
    for(j=0;j<NMOL;j++){
      for (xi=0;xi<3;xi++) rc[xi]=0.0; //initialization
      for (k=0;k<NAT;k++) { //Center of mass loop
	a=fscanf(fcnf, "%*5d%*5s%*5s%*5d%f%f%f%*f%*f%*f\n",rr+0,rr+1,rr+2); //only interested in x,y,z
        if(ferror(fcnf) ) printf("ferror found!\n");
        if(feof(fcnf) ) printf("feof found!\n");
        if(a!=3) {
          fprintf(stderr,"Error in reading molecule %d, atom %d from %s:\n fscanf returned %d, read %f %f %f\n",j+1,k+1,cnffile,a,rr[0],rr[1],rr[2]);
          fprintf(stderr," Previously read position:\t %f %f %f\n",foo[0],foo[1],foo[2]);
          fprintf(stderr,"You should check file %s around line %d\n\n",cnffile, 2+(j+1)*NAT+(k+1));
          exit(1);
        }
        // Molecule's center of mass with PBC
	for (xi=0;xi<3;xi++) {
          foo[xi]=rr[xi]; //backup for error
          if(k==0) r1[xi]=rr[xi]; // fix atom 1 position
          d = rr[xi]-r1[xi];
          d-=L[xi]*floor(d/L[xi] + 0.5); // apply PBC to molecules that cross the boundary at a fixed frame
          rc[xi]+=(double)(r1[xi] + d);
        }
      }
      // 3) LOOP ON x,y,z
      for (xi=0;xi<3;xi++) {
	idx=NMOL*3*i+3*j+xi;
	r[idx]=rc[xi]/NAT;
	if(i>0){
	  idx_old=idx-(NMOL*3);
	  dr=r[idx]-r[idx_old];     //displacement in t1-first_t1, w.o. PBC
	  //if(idx==0) fprintf(stdout,"%f %f ",t,Dr[idx]+dr);
	  dr-=L[xi]*floor(dr/L[xi] + 0.5);  //PBC
	  r[idx]=r[idx_old]+dr;     //virtual position at t1, with PBC wrt first_t1
	  // 4) LOOP ON t0<t1 (invalid for the first t1)
	  for(a=0;a<i;a++) {
	    dt=t1-ts[i0+a*NPC];
	    idx_old=NMOL*3*a+3*j+xi;
	    dr = r[idx]-r[idx_old];
	    dr2=dr*dr;
	    dr4=dr2*dr2;
	    ti = floor(dt/tcycle);
	    msd[ti]+=dr2;
	    msd2[ti]+=dr4;
	    if((xi%3)==0) Navg[ti]++; //don't divide <r^2> by 3
	  } // 4)
	} // 3)
      } // 2)
    } // 1)
    fclose(fcnf);
    fprintf(stderr," t1=%10lf (%d/%d); done from %s\n",t1,i+1,ncycles,cnffile);
  }
  /*
  fprintf(stderr,"\nAveraging (r(t)-r(t0))^2 over t0,particles and coordinates...");
  for(i=0;i<ncycles;i++) {
    t0=ts[i0+i*NPC];
    for(a=0;a<=ncycles-i;a++) {
      dt=ts[i0+(a+i)*NPC]-t0;
      for(j=0;j<NMOL;j++) {
	for(xi=0;xi<3;xi++) {
	  idx_old = NMOL*3*i+3*j+xi;
	  idx = NMOL*3*(a+i)+3*j+xi;
	  dr = r[idx]-r[idx_old];
	  dr2=dr*dr;
	  dr4=dr2*dr2;
	  ti = floor(dt/tcycle +0.5);
	  msd[ti]+=dr2;
	  msd2[ti]+=dr4;
	  Navg[ti]++;
	}
      }
    }
  }
  */
  sprintf(outfile,"msd_nat%d-lin.dat",NAT);
  fprintf(stderr," done.\n\n Saving MSD on %s:\n",outfile);
  fout=fopen(outfile,"w");
  for(i=1;i<ncycles;i++) {
    dt=ts[i0+i*NPC]-ts[i0];
    ti=floor(dt/tcycle);
    msd[ti]/=Navg[ti];
    msd2[ti]=sqrt( msd2[ti]/Navg[ti]-(msd[ti]*msd[ti]) );
    if(Navg[ti]>NMOL) msd2[ti]/=sqrt(Navg[ti]/NMOL-1); //average over cycles, not atoms
    fprintf(fout, "%18le %18le %18le\n",dt,msd[ti],msd2[ti]); // sigma/total n.of pts of the average
    fprintf(stderr," dt=%10f; MSD(dt)=%10f +/- %10f (averaged over %3d t0's)\n",dt,msd[ti],msd2[ti],Navg[ti]/NMOL);
  }
  fclose(fout);

  /**************** CENTER OF MASS M.S.D *******************/

  fprintf(stderr," done\n\nEvaluating Center of Mass trajectory... ");
  sprintf(outfile,"trjcm_nat%d-lin.dat",NAT);
  fout=fopen(outfile,"w");
  for(i=0;i<ncycles;i++) {
    fprintf(fout,"%lf",ts[i0+i*NPC]);
    for(xi=0;xi<3;xi++) {
      idx_old=NMOL*3*i+3*0+xi; //store center of mass in particle idx_old (j==0)
      for(j=1;j<NMOL;j++) { //add to it other r's
	idx = NMOL*3*i+3*j+xi;
	r[idx_old] += r[idx];
      }
      r[idx_old]/=NMOL;
      fprintf(fout," %lf",r[idx_old]);
    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  fprintf(stderr,"done. Saved on %s\n",outfile);

  fprintf(stderr,"\nAveraging center of mass (r(t0+dt)-r(t0))^2 over t0...");
  for(ti=0;ti<maxncycles;ti++) {
    msd[ti]=0.0;
    msd2[ti]=0.0;
    Navg[ti]=0;
  }
  for(i=0;i<ncycles;i++) {
    t1=ts[i0+i*NPC];
    for(xi=0;xi<3;xi++) {
      idx = NMOL*3*i+3*0+xi;
      if(i>0) {
	for(a=0;a<i;a++) {
	  dt=t1-ts[i0+a*NPC];
	  idx_old = NMOL*3*a+3*0+xi;
	  dr = r[idx]-r[idx_old];
	  dr2=dr*dr;
	  dr4=dr2*dr2;
	  ti = floor(dt/tcycle);
	  if(ti>=maxncycles) fprintf(stderr,"ERROR: i=%d dt=%lf ti=%d exceeds maxncycles\n",i,dt,ti);
	  msd[ti]+=dr2;
	  msd2[ti]+=dr4;
	  if((xi%3)==0) Navg[ti]++; //don't divide by 3 <r^2(t)>
	}
      }
    }
  }
  sprintf(outfile,"msdcm_nat%d-lin.dat",NAT);
  fout=fopen(outfile,"w");
  fprintf(stderr," done.\n\n Saving MSDcm on %s\n",outfile);
  for(i=1;i<ncycles;i++) {
    dt=ts[i0+i*NPC]-ts[i0];
    ti=floor(dt/tcycle);
    if(ti>=maxncycles) fprintf(stderr,"ERROR: i=%d dt=%lf ti=%d exceedes maxncycles\n",i,dt,ti);
    msd[ti]/=Navg[ti];
    msd2[ti]/=Navg[ti];
    msd2[ti]=sqrt(msd2[ti]-(msd[ti]*msd[ti]));
    if(Navg[ti]>1) msd2[ti]/=sqrt(Navg[ti]-1);
    fprintf(fout, "%18le %18le %18le\n",dt,msd[ti],msd2[ti]); // sigma/total n.of pts of the average
    fprintf(stderr," dt=%10f; MSDcm(dt)=%10f +/- %10f (averaged over %3d t0's)\n",dt,msd[ti],msd2[ti],Navg[ti]);
  }
  fclose(fout);
  fprintf(stderr,"done\n");

  free(r);
  free(msd);
  free(msd2);
  free(Navg);
  return 0;
}

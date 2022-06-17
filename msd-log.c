#include "lib.h"

#define JREF 0    // index of the molecule whose dr^2 will be live-printed
#define XIREF 0   // index of the coordinate x,y,z which will be live-printed


void printErr(char *myname) {
  fprintf(stderr, "\nUsage: %s [t0] [t1] \n\nQuesto programma calcola il MSD(t) a partire dai frame corrispondenti contenuti nel file lista.dat.\n t0\ttempo iniziale (opzionale), in ps\n t1\ttempo finale (opzionale), in ps\n\n", myname);
  exit(1);
}

int main(int argc, char **argv) {
  int a,i,i0,i1,j,k,xi,nt,Ncnf,idx, steps[MAXCNF], ncycles;
  float r[3], r1[3], d, foo[3];
  double t,t0,t1,ts[MAXCNF],tcycle, msd[NPC],msd2[NPC],dr,Dr2,rc[3], *r0, r2,r4, msdcm[NPC],msdcm2[NPC], DRcm;
  char str[10], cnffile[50], outfile[30];
  FILE *fcnf,*fout, *foutcm;
  Ncnf = getTimes(steps, ts); //returns n. of configurations, their time (ps), and initializes cnfroot, N, L !!!!!!!!!
  i0=0;
  i1=Ncnf-1;
  fprintf(stderr,"Available time window: [%lf, %lf] ps\n",ts[i0],ts[i1]);
  if(argc>3) printErr(argv[0]);
  else if (argc>=2) {
    t=atof(argv[1]);
    while(ts[i0]<t && i0<Ncnf-1) i0++;
    i0=(i0/NPC)*NPC; //prendi il tempo zero piu' vicino da sotto
    if (argc==3) {
      t=atof(argv[2]);
      while(ts[i1]>t && i1>0) i1--;
      i1=((i1-i0)/NPC)*NPC+i0;
    }
  }
  nt=i1-i0+1;
  t0=ts[i0];
  t1=ts[i1];
  if(Ncnf<NPC) {
    fprintf(stderr,"Error: available data is less than 1 log-cycle (>%f ps)\n",ts[Ncnf-1]-ts[0]);
    exit(1);
  }
  tcycle = ts[i0+NPC-1]-t0;
  ncycles=nt/NPC;
  fprintf(stderr,"You chose: \t\t[%lf, %lf] ps\n => %d cycles (%f ps each)\n",t0,t1, ncycles,tcycle);
  if(ncycles==0) {
    fprintf(stderr,"Error: selected time is less than 1 log-cycle.\n");
    exit(1);
  }

  r0=(double *)malloc(NMOL*3*sizeof(double));
  for(i=0;i<NPC;i++) msd[i]=msd2[i]=0.0;
//  fprintf(stderr, "\nL=%f nm\n\nCoordinate %d of molecule idx.%d will be live-printed\n",L,JREF,XIREF);

  fprintf(stderr, "Simultaneous computation of MSD and MSDcm started\n\n");
  // 1) LOOP ON t0
  for(a=0;a<ncycles;a++) {
    t0=ts[i0+a*NPC];
    // 2) LOOP ON t1 -> delta t
    for(i=0;i<NPC;i++) {
      r2=r4=0.0;
      t=ts[i0+a*NPC+i]-t0;
      sprintf(str,"-%d",steps[i0+a*NPC+i]);
      strcpy(cnffile,cnfroot);
      strcat(cnffile,str);
      fcnf=fopen(cnffile,"r");
      read_N_L(fcnf, &N, L);
      // 3) LOOP ON MOLECULES
      DRcm=0.0;
      for(j=0;j<NMOL;j++){
	for (xi=0;xi<3;xi++) rc[xi]=0.0; //initialization
	for (k=0;k<NAT;k++) { //Center of mass loop
	  xi=fscanf(fcnf, "%*5d%*5s%*5s%*5d%f%f%f%*f%*f%*f\n",r+0,r+1,r+2); //only interested in x,y,z
          if(ferror(fcnf) ) printf("ferror found!\n");
          if(feof(fcnf) ) printf("feof found!\n");
          if(xi!=3) {
            fprintf(stderr,"Error in reading molecule %d, atom %d from %s:\n fscanf returned %d, read %f %f %f\n",j+1,k+1,cnffile,xi,r[0],r[1],r[2]);
            fprintf(stderr," Previously read position:\t %f %f %f\n",foo[0],foo[1],foo[2]);
            fprintf(stderr,"You should check file %s around line %d\n\n",cnffile, 2+(j+1)*NAT+(k+1));
            exit(1);
          }
          // Molecule's center of mass with PBC
	  for (xi=0;xi<3;xi++) {
            foo[xi]=r[xi];
            if(k==0) r1[xi]=r[xi]; // fix atom 1 position
            d = r[xi]-r1[xi];
            d-=L[xi]*floor(d/L[xi] + 0.5); // apply PBC to molecules that cross the boundary at a fixed frame
            rc[xi]+=(double)(r1[xi] + d);
          }
	}
	// 4) LOOP ON x,y,z
	for (xi=0;xi<3;xi++) {
	  rc[xi]/=NAT;
	  idx=3*j+xi;
	  if(i==0) {
	    r0[idx]=rc[xi];//if Dt is 0, just save rc into r0
//	    if(j==JREF && xi==XIREF) fprintf(stderr, "\nloaded position %f \t from %s\n",r0[idx],cnffile);
	  }
	  else {
	    dr=rc[xi]-r0[idx];    //displacment from t=0 (it's small enough not to make 2 round trips around the box)
	    //if(idx==999) fprintf(stdout,"%f %f ",t,Dr[idx]+dr);
	    dr-=L[xi]*floor(dr/L[xi] + 0.5); //apply PBC to c.o.m. displacement in time
            DRcm+=dr;
	    Dr2=dr*dr;    //average Dr^2 over particles
//	    if(j==JREF && xi==XIREF) fprintf(stderr,"dt=%10f; dr^2=%10g \t from %s\n",t,Dr2,cnffile);
	    r2+=Dr2;
	    r4+=(Dr2*Dr2);
	  }
	} // 4)
      } // 3) end of molecules
      fclose(fcnf);
      msd[i]+=(r2/NMOL);
      msd2[i]+=(r4/NMOL);
      DRcm/=NMOL;
      msdcm[i]+=(DRcm*DRcm);
      msdcm2[i]+=(DRcm*DRcm)*(DRcm*DRcm);
      //fprintf(stderr,"\rt0=%10lf (%d/%d); t=%10lf (%d/%d) done from %s",t0,a+1,ncycles,t,i+1,NPC,cnffile);
    } // 2) end of t1 at fixed t0
    fprintf(stderr,"\nt0=%10lf (%d/%d) done\n",t0,a+1,ncycles);
  } //1) end of t0
  
  fprintf(stderr,"\n");
  sprintf(outfile,"msd_nat%d-log.dat",NAT);
  fout=fopen(outfile,"w");
  sprintf(outfile,"msdcm_nat%d-log.dat",NAT);
  foutcm=fopen(outfile,"w");
  fprintf(stderr, "\nSaving MSD on msd_natX-log.dat and MSDcm on msdcm_natX-log.dat\n");
  for(i=0;i<NPC;i++) {
    t=ts[i0+i]-ts[i0];
    msd[i]/=ncycles;
    msd2[i]=sqrt( msd2[i]/ncycles-(msd[i]*msd[i]) );
    msdcm[i]/=ncycles;
    msdcm2[i]=sqrt( msdcm2[i]/ncycles-(msdcm[i]*msdcm[i]) );
    if(ncycles>1) msd2[i]/= sqrt((ncycles-1));
    if(i>0){ //don't print dt=0
      fprintf(fout, "%18le %18le %18le\n",t,msd[i],msd2[i]);
      fprintf(foutcm, "%18le %18le %18le\n",t,msdcm[i],msdcm2[i]);
    }
    fprintf(stderr,"(%d/%d) dt=%10f; MSD(dt)=%10f +/- %10f; MSDcm(dt)=%10f +/- %10f\n", i,NPC-1,t,msd[i],msd2[i],msdcm[i],msdcm2[i]);
  }
  fclose(fout);
  fprintf(stderr,"\n\n");
  //free(r0);
  return 0;
}

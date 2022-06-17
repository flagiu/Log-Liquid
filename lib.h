#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<complex.h>

#define NAT 1          // atoms per molecule
int N;
//#define N 708588       // n. of atoms
float L[3];  // box size (nm), read from configuration
//#define L 43.32963     // box size (nm)
#define MAXNMAX 200    // nmax over which the exp array is too big (experimented on my Acer Aspire E15, 4GB RAM)
#define MAXCNF 5000    // max n. of configuration files
#define MAXQARR 300    // max n. of qvectors with same q
#define DN 0.5         // chosen mesh size in QVECTORSxxx/ files
#ifndef MESHFACT
#define MESHFACT 1     //(positive integer) increases the mesh size
#endif
#ifndef NPC
#define NPC 50         // log-spaced times' cycle length
#endif
#define LOG 0          // timing mode: logarithmic or linear
#define LIN 1

#define NMOL (N/NAT)   //number of molecules
#define sqrtNMOL (sqrt((double)(NMOL)))

double Dq[3];		//minimum resolution in qx,qy,qz due to box periodicity
double dq[3];		//minimum resolution for q chosen in QVECTORSxxx/
//#define QNMAX 5         // maximum q-vector modulus (in units of Dq)

unsigned int qnmax; // maximum absolute wavenumber to be calculated along any direction (it will be defined in isf.c/main), in units of Dq!!
char cnfroot[30];   // root name for configuration files, usually ../CnfLT???

#define N1 NMOL
#define N2 (qnmax+1) // range of wavenumbers to be calculated (0,1,..,qnmax)
#define N3 3 //3 dimensions
#define N4 2 //real and imag
#define EXPN (N1*N2*N3*N4)
#define EXPidx(i,qn,xi,j) (i*N2*N3*N4 + qn*N3*N4 + xi*N4 + j)
#define EXPi(idx) (idx%(N2*N3*N4))
#define EXPqn(idx,i) ((idx-i*(N2*N3*N4))%(N3*N4))
#define EXPx(idx,i,qn) ((idx-i*N2*N3*N4-qn*N3*N4)%N4)
#define EXPj(idx,i,qn,xi) (idx-i*N2*N3*N4-qn*N3*N4-xi*N4)

// Real and Imag part of (a0+i*b0)*(a1+i*b1)
//#define RE_2PRODUCT(a0,b0,a1,b1) (a0*a1 - b0*b1)
//#define IM_2PRODUCT(a0,b0,a1,b1) (a0*b1 + b0*a1)
// Real and Imag part of (a0+i*b0)*(a1+i*b1)*(a2+i*b2)
//#define RE_3PRODUCT(a0,b0,a1,b1,a2,b2) (a0*a1*a2 - a0*b1*b2 - a1*b2*b0 - a2*b0*b1)
//#define IM_3PRODUCT(a0,b0,a1,b1,a2,b2) (-b0*b1*b2 + b0*a1*a2 + b1*a2*a0 + b2*a0*a1)

// read N and L and set the pointer to be ready to read coordinates
void read_N_L(FILE* fp, int *Nptr, float L[3]) {
	char c;
	int res;
	fseek(fp, -1, SEEK_END);
	do{ // find the last '\n'
		fseek(fp, -2, SEEK_CUR);
		c = fgetc(fp);
	} while (c!='\n');
	fseek(fp, 1, SEEK_CUR); // go one step ahead to find the last line
	res=fscanf(fp, "%f%f%f\n", L,L+1,L+2); // scan Lx,Ly,Lz
	fseek(fp, 0, SEEK_SET); //go back to beginning
	do{ // find the first '\n' (note: you don't need fseek to go ahead
		c = fgetc(fp);
	} while (c!='\n');
	res=fscanf(fp, "%d\n", Nptr); // scan the number of particles
//	fprintf(stderr, "\nN = %d atoms\nL[] = %f %f %f nm\n\n",*Nptr, L[0],L[1],L[2]);
}

// store coordinates (r) or exp(i*q*r), depending on onfly, from a .gro file into the ptr muldimensional array
void get_gro(double *ptr, int step, int qnmax) {
  FILE *fcnf;
  char cnffile[40]={0},str[10];
  int i,qn,xi,j, idx;
  float r[3], r1[3], d, foo[3];
  double rc[3], arg;
  sprintf(str,"-%d",step);
  strcat(cnffile,cnfroot);
  strcat(cnffile,str);//sprintf(cnffile,"../Cnfem95-%d",step);
  if(qnmax==-1) fprintf(stderr, "Storing coordinates from %s ...",cnffile);
  else fprintf(stderr, "Storing exp-array from %s ...",cnffile);
  fcnf=fopen(cnffile,"r");
  read_N_L(fcnf, &N, L);
  for(xi=0;xi<3;xi++) {
    Dq[xi] = 2*M_PI/L[xi];  //minimum resolution for qx,qy,qz (physical constraints)
    dq[xi] = DN*Dq[xi];     //minimum " " " chosen in QVECTORSxxx/
  }
  // 1) MOLECULES LOOP
  for(i=0;i<NMOL;i++) {
    for (xi=0;xi<3;xi++) rc[xi]=0.0; //initialization
    // 1.1) Center of mass calculation loop
    for (j=0;j<NAT;j++) {
      xi=fscanf(fcnf, "%*5d%*5s%*5s%*5d%f%f%f%*f%*f%*f\n",r+0,r+1,r+2); //only interested in x,y,z
      if(ferror(fcnf) ) printf("ferror found!\n");
      if(feof(fcnf) ) printf("feof found!\n");
      if(xi!=3) {
	fprintf(stderr,"Error in reading molecule %d, atom %d from %s:\n fscanf returned %d, read %f %f %f\n",i+1,j+1,cnffile,xi,r[0],r[1],r[2]);
	fprintf(stderr," Previously read position:\t %f %f %f\n",foo[0],foo[1],foo[2]);
	fprintf(stderr,"You should check file %s around line %d\n\n",cnffile, 2+(i+1)*NAT+(j+1));
	exit(1);
      }
      // Molecule's center of mass with PBC
      for (xi=0;xi<3;xi++) {
	foo[xi]=r[xi];
	if(j==0) r1[xi]=r[xi]; // fix atom 1 position
	d = r[xi]-r1[xi];
	d-=L[xi]*floor(d/L[xi] + 0.5); // apply PBC to molecules that cross the boundary at a fixed frame
	rc[xi]+=(double)(r1[xi] + d);
      }
    }
    // 1.2) COORDINATES LOOP
    for(xi=0;xi<3;xi++) {
      rc[xi]/=NAT; //normalize center of mass coord.
      // 1.2.1) if qnmax is not given (aka -1 is passed), STORE only rc
      if(qnmax==-1) ptr[i*3+xi]=rc[xi];
      else {
	// 1.2.2) QVECTOR LOOP: if qnmax is given, STORE every exp(i*q*r)
	for (qn=0;qn<=qnmax;qn++) {
	  idx=EXPidx(i,qn,xi,0);
	  arg=Dq[xi]*qn*rc[xi]; //remember, qnmax is in units of Dq
	  ptr[idx]  =cos(arg);
	  ptr[idx+1]=sin(arg);
	  //fprintf(stdout,"%d %18le %18le\n",idx,exp[idx],exp[idx+1]);
	} // 1.2.2)
      }
    } //1.2)
    //if (((i+1)%100000)==0) fprintf(stderr, "\r %d/%d molecules evaluated.",i+1,NMOL);
  } // 1) 
  fprintf(stderr," done.");
  fclose(fcnf);
}
void get_exp(double *exp, int step, int qnmax) {
  get_gro(exp,step,qnmax);
}
void get_coords(double *r, int step) {
  get_gro(r,step,-1);
}
void calcRHOfromE(double rho[2], double *e, int qvec[3], int printexp) {
  int i,xi, idx;
  double c[3],s[3], Re=0.0,Im=0.0; //cos,sin and Re,Im parts of rho(q,t)
  for(i=0;i<NMOL;i++) {
    for(xi=0;xi<3;xi++) {
      idx=EXPidx(i,abs(qvec[xi]),xi,0);
      c[xi]=e[idx];
      s[xi]=( qvec[xi]<0 ? -e[idx+1] : e[idx+1] );
      //debugging
      if (printexp && qvec[xi]!=0) fprintf(stdout,"%d %18le %18le\n",i+1,c[xi],s[xi]);
    }
    Re+=(  c[0]*c[1]*c[2] - s[0]*s[1]*c[2] - s[1]*s[2]*c[0] - s[2]*s[0]*c[1] );
    Im+=( -s[0]*s[1]*s[2] + s[0]*c[1]*c[2] + s[1]*c[2]*c[0] + s[2]*c[0]*c[1] );
  }
  rho[0]=Re/sqrtNMOL;
  rho[1]=Im/sqrtNMOL;
}
void calcRHOfromR(double rho[2], double *r, int qvec[3], int printexp) {
  int i,xi;
  double arg, c[3],s[3], Re=0.0,Im=0.0; //cos,sin and Re,Im parts of rho(q,t)
  for(i=0;i<NMOL;i++) {
    for(xi=0;xi<3;xi++) {
      arg=Dq[xi]*qvec[xi]*r[i*3+xi];
      c[xi]=cos(arg);
      s[xi]=sin(arg);
      //debugging
      if (printexp && qvec[xi]!=0) fprintf(stdout,"%d %18le %18le\n",i+1,c[xi],s[xi]);
    }
    Re+=(  c[0]*c[1]*c[2] - s[0]*s[1]*c[2] - s[1]*s[2]*c[0] - s[2]*s[0]*c[1] );
    Im+=( -s[0]*s[1]*s[2] + s[0]*c[1]*c[2] + s[1]*c[2]*c[0] + s[2]*c[0]*c[1] );
  }
  rho[0]=Re/sqrtNMOL;
  rho[1]=Im/sqrtNMOL;
}
void calcRHORHO(double fqt[2], double *e0, double *e1, int qvec[3], double rho0[2], double rhot[2], int printexp, int onfly) {
  if(onfly) {
    calcRHOfromR(rho0,e0,qvec, 0);
    calcRHOfromR(rhot,e1,qvec, printexp);
  }
  else {
    calcRHOfromE(rho0,e0,qvec, 0);
    calcRHOfromE(rhot,e1,qvec, printexp);
  }
  // Re(rho(t)rho(0)*) = Re(t)*Re(0) + Im(t)*Im(0)
  fqt[0]=rhot[0]*rho0[0]+rhot[1]*rho0[1];
  // Im(rho(t)rho(0)*) = Re(t)*Im(0) - Im(t)*Re(0)
  fqt[1]=rhot[0]*rho0[1]-rhot[1]*rho0[0]; 
}

void calcFself(double Fs[2], double *e0, double *e1, int qvec[3], int onfly) {
  int i,idx,j;
  double c,s;
  double _Complex exp, fs = 0.0 + 0.0*I; // without complex type here we would go crazy...
  if(onfly) {
    for(i=0;i<NMOL;i++) {
      exp = 1.0 + 0.0*I;
      for(j=0;j<3;j++) {
	exp *= cexp(-1*I * Dq[j]*qvec[j]*e0[i*3+j]); // exp(-q*r_i) at t=t0
        exp *= cexp(   I * Dq[j]*qvec[j]*e1[i*3+j]); // exp(q*r_i) at t=t1
      }
      fs += exp;
    } 
  }
  else {
    for(i=0;i<NMOL;i++) {
      exp = 1.0 + 0.0*I;
      for(j=0;j<3;j++) {
	idx=EXPidx(i,abs(qvec[0]),j,0);
	c=e0[idx];
	s=( qvec[0]<0 ? -e0[idx+1] : e0[idx+1] );
	exp *= (c - s*I); // exp(-q*r_i) at t=t0
	c=e1[idx];
	s=( qvec[0]<0 ? -e1[idx+1] : e1[idx+1] );
	exp *= (c + s*I); // exp(q*r_i) at t=t1
      }
      fs += exp;
    }
  }
  Fs[0] = creal(fs)/NMOL;
  Fs[1] = cimag(fs)/NMOL;
}


//////////////////////////////////////////////////////////////////////////////

// reads configurations in lista.dat, deduces their root name, saves the corresponding run.sh step and time (ps), INITIALIZES N and L, and returns the number of files
int getTimes(int *steps, double *ts) {
  int i, read=0;
  char cnfname[30];
  FILE *fp;
  fp=fopen("lista.dat","r");
  for(i=0;read!=-1;i++) {
    read=fscanf(fp, "%s %d %lf\n",cnfname,steps+i,ts+i);
  }
  fclose(fp);
  fprintf(stderr, "\nNAT = %d grouped atoms\nNPC = %d points per time cycle\nInitialized N,L,cnfroot from last configuration %s\n",NAT,NPC,cnfname);
  // initialize N and Lx,Ly,Lz !!!
  fp=fopen(cnfname,"r");
  read_N_L(fp, &N, L);
  fclose(fp);
  // take the root of the (last) name
  strcpy(cnfroot, strtok(cnfname,"-") );
  fprintf(stderr, " N = %d atoms\n L = %f %f %f nm\n cnfroot = %s\n", N,L[0],L[1],L[2],cnfroot);
  return i-1; //return number of lines
}

// for checkpoint restart; we save t0, Dt and wavenumber
void restart(char *logfile, int *t0, int *Dt, int *qn) {
  FILE *fp;
  int foo;
  if ( (fp=fopen(logfile,"r")) != NULL) {
    foo=fscanf(fp,"%d\n",t0);
    foo=fscanf(fp,"%d\n",Dt);
    foo=fscanf(fp,"%d\n",qn);
    fclose(fp);
    fprintf(stderr, "Restarting from %s: ", logfile);
  } else {
    fprintf(stderr, "%s not found: ", logfile);
    *t0=*qn=0;
    if(Dt!=NULL) *Dt=0; //it is null when called from struct.c
  }
  fprintf(stderr,"\t t0=%d, ",*t0);
  if(Dt != NULL)  fprintf(stderr,"Dt=%d, ",*Dt);
  fprintf(stderr,"nqvec=%d\n",*qn);
}

// checkpoint saving into a .logc file
void savelog(char *logfile, int t0, int Dt, int qn) {
  FILE *fp;
  fp=fopen(logfile,"w");
  fprintf(fp,"%d\n",t0);
  fprintf(fp,"%d\n",Dt);
  fprintf(fp,"%d\n",qn);
  fflush(fp);
  fclose(fp);
}

// reads the list of qvectors [nx,ny,nz] from QVECTORxxx/qvector.xxx into qarr[][]
int getQvectors(char *qfile, int qarr[MAXQARR][3]) {
  FILE *fp;
  int i,read=0;
  fp=fopen(qfile, "r");
  if(qarr==NULL)
    for(i=0;read!=EOF;i++)
      read=fscanf(fp,"%*d %*d %*d\n");//if NULL argument: only count lines
  else
    for(i=0;read!=EOF;i++)
      read=fscanf(fp, "%d %d %d\n",*(qarr+i)+0,*(qarr+i)+1,*(qarr+i)+2);//qarr[i]+0,qarr[i]+1,qarr[i]+2);
  fclose(fp);
  return i-1;
}

// reads the min and max wavenumber in QVECTORxxx/qvector.xxx (for debugging)
void getMinMaxWavenumber(char *qfile, unsigned int out[2]) {
  FILE *fp;
  int i,j,read=0, n[3],m;
  unsigned int min=10000, max=0;
  fp=fopen(qfile, "r");
  for(i=0;read!=EOF;i++) {
    read=fscanf(fp, "%d %d %d\n",n,n+1,n+2);
    for(j=0;j<3;j++) {
      m=abs(n[j]); //in absolute value!
      if(m<min) min=m;
      if(m>max) max=m;
    }
  }
  fclose(fp);
  out[0]=min;
  out[1]=max;
  return;
}

// live printing of the status; it works better on a wide window terminal
void updateTerminal(double t0, int i0, int N0, double Dt, int iDt, int NDt, int qvec[3], int iq, int Nq) {
  fprintf(stderr, "\r t0=%f (%d/%d); Dt=%f (%d/%d); qvector=[%d,%d,%d] (%d/%d); done.\033[K", t0,i0,N0,Dt,iDt,NDt,qvec[0],qvec[1],qvec[2],iq,Nq);
  if (iq==Nq) fprintf(stderr,"\n");
}

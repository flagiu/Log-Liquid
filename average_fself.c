#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.0001
#define MAXLINES 100000000
//!!! the input file MUST be ordered by the 2nd or 3rd column (i.e. time interval) !!!

void printaverage(double dt, double dt2, double R, double R2, double I, double I2, int N) {
  double mt,st,mR,mI,sR,sI, M=(double)N;
  mt=dt/M;
  st=sqrt( (dt2/M - mt*mt)/M );
  mR=R/M;
  sR=sqrt( (R2/M - mR*mR)/M );
  mI=I/M;
  sI=sqrt( (I2/M - mI*mI)/M );
  printf("%18le %18le %18le %18le %18le %18le %d\n",mt,st, mR,sR, mI,sI, N);
  //fprintf(stderr,"\nstatistics for dt=%10lf done over %d pts\n",mt,N);
}

void exitErr(char* myname) {
  fprintf(stderr,"ERROR: wrong input.\n\nUsage: %s <input file> [n. of lines]\n",myname);
  exit(1);
}

int main(int argc, char **argv) {
  int j=0,newj,foo, firstrow=1,N, res, count=0, Nlines=MAXLINES;
  double t,x,y,dt,dt2,R,R2,I,I2;
  FILE *fp;
  if (argc!=2 && argc!=3) exitErr(argv[0]);
  if (argc==3) {
    Nlines=atoi(argv[2]);
    if (Nlines<=0) exitErr(argv[0]);
  }
  dt=dt2=R=I=R2=I2=0.0;
  N=0; //number of pts for the average
  fp=fopen(argv[1], "r");
  //fprintf(stderr,"\n");
  if (fp != NULL) {
    do {
      res = fscanf(fp, "%d %d %18le %d %d %d %18le %18le\n",
		   &foo,&newj,&t,&foo,&foo,&foo,&x,&y);
      if (firstrow) firstrow=0;
      else if (newj!=j) {
	//calculate and print average for the block that just ended
	printaverage(dt,dt2,R,R2,I,I2,N);
	// go ahead with the next dt-block
	//fprintf(stderr,"\nj=%10d \t N=%10d \t newj=%10d\n",j,N,newj); 
	j=newj;
	dt=dt2=R=I=R2=I2=0.0;
	N=0;
      }
      dt+=t;
      dt2+=(t*t);
      R+=x;
      R2+=(x*x);
      I+=y;
      I2+=(y*y);
      N++;
      count++;
      //fprintf(stderr,"\rcount=%8d dt-index=%5d res=%2d",count,newj,res);
    } while (res!=EOF && count<=Nlines+1);
    //last block
    printaverage(dt,dt2,R,R2,I,I2,N);
    fclose(fp);
    //fprintf(stderr,"\n");
    if (count==Nlines+2) fprintf(stderr,"ERROR: maximum number of lines exceeded (%d/%d) without reaching EOF.\n\n",count-1,Nlines);
  } else {
    fprintf(stderr, "ERROR: could not open %s\n",argv[1]);
    exit(1);
  }
  return 0;
}

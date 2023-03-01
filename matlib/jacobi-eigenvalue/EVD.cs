using System;
using static System.Math;
using static System.Console;
public class EVD{

public vector e;
public matrix V=null;
public int sweeps;

public static void timesJ(matrix A, int p, int q, double c, double s){
        for(int i=0;i<A.size1;i++){
                double new_aip=c*matrix.get(A,i,p)-s*matrix.get(A,i,q);
                double new_aiq=s*matrix.get(A,i,p)+c*matrix.get(A,i,q);
                matrix.set(A,i,p,new_aip);
                matrix.set(A,i,q,new_aiq);
                }

	}

public static void Jtimes(matrix A, int p, int q, double c,double s){
        for(int j=0;j<A.size2;j++){
                double new_apj= c*matrix.get(A,p,j)+s*matrix.get(A,q,j);
                double new_aqj=-s*matrix.get(A,p,j)+c*matrix.get(A,q,j);
                matrix.set(A,p,j,new_apj);
                matrix.set(A,q,j,new_aqj);
                }
	}//Jtimes

public EVD(matrix A, bool withV=true){// A is destroyed
	int n=A.size1;
	if(withV) V = matrix.id(n);
	e = new vector(n);
	bool changed;
	sweeps=0;
	do{
		changed=false;
		sweeps++;
		for(int p=0;p<n-1;p++)
		for(int q=p+1;q<n;q++) {
			double apq=matrix.get(A,p,q);
			double app=matrix.get(A,p,p);
			double aqq=matrix.get(A,q,q);
			double theta=0.5*Atan2(2*apq,aqq-app);
			double c=Cos(theta),s=Sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			if(new_app!=app || new_aqq!=aqq) {
				changed=true;
				timesJ(A,p,q,c,s);
				Jtimes(A,p,q,c,-s);
				if(withV) timesJ(V,p,q,c,s);
			}
		}
	}while(changed);
	for(int i=0;i<n;i++)e[i]=A[i,i];
}//ctor


}//EVD

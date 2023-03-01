using System;
using static System.Console;
using static System.Math;
class main {

public static double data(double x){
//	if(x<0)return 0; else return 1;
	return 1/(1+x*x);
	}

public static double polinterp(double[] x,double[] y,double z){
	int n=x.Length;
	double s=0,p;
	for(int i=0;i<n;i++) {
		p=1;
		for(int k=0;k<n;k++) if(k!=i) p*=(z-x[k])/(x[i]-x[k]);
		s+=y[i]*p; }
	return s;
	}	

public static void Main(){

int i, n=9, N=200;
double[] x = new double[n];
double[] y = new double[n];

double z;
double a=-6, b=6;
for(i=0; i<n; i++){
	x[i]=a+(b-a)*i/(n-1);
	y[i]=data(x[i]);
	WriteLine($"{x[i]} {y[i]}");
	}
Write("\n\n");

double step=(x[n-1]-x[0])/(N-1);
for (z=x[0], i=0; i<N; z=x[0]+(++i)*step){
	WriteLine($"{z} {polinterp(x,y,z)} {matlib.berrut(x,y,z)}");
	}

}//Main

}//main

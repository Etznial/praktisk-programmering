using System;
using static System.Console;
using static System.Math;
class main {

static void Main(){

int n=5, N=200;
double[] x = new double[n];
double[] y = new double[n];

int i;
for (i=0; i<n; i++){
	x[i]=2*PI*i/(n-1);
	y[i]=Sin(x[i]);
	WriteLine("{0:g6} {1:g6}",x[i],y[i]);
	}

Write("\n\n");

var cs = new cspline(x,y);
double z, step=(x[n-1]-x[0])/(N-1);
for (z=x[0], i=0; i<N; z=x[0]+(++i)*step){
	WriteLine($"{z} {Sin(z)} {cs.eval(z)}");
	}
Write("\n\n");
for (z=x[0], i=0; i<N; z=x[0]+(++i)*step){
	WriteLine($"{z} {Cos(z)} {cs.deriv(z)}");
	}
Write("\n\n");
for (z=x[0], i=0; i<N; z=x[0]+(++i)*step){
	WriteLine($"{z} {1-Cos(z)} {cs.integ(z)}");
	}

}//Main

}//main

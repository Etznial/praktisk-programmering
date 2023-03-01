using System;
using static System.Console;
using static System.Math;
using System.Collections;
using System.Collections.Generic;

public static partial class matlib{

static bool approx(double a,double b, double acc=1e-6, double eps=1e-6){
	if( Abs(a-b)<acc )return true;
	if( Abs(a-b)<eps*Max(Abs(a),Abs(b)) ) return true;
	return false;
}

static Func<double,vector,vector>
F=delegate(double x, vector y){
	return new vector(Exp(-x*x)*2/Sqrt(PI));
		};

public static double erf(double z){
	if(z<0)return -erf(-z);
	double a=0;
	vector ya=new vector(0.0);
	double h=0.05,acc=1e-3,eps=1e-3;
	vector y = ODE.driver(F,a,ya,z,acc:acc,eps:eps,h:h);
	return y[0];
	}
}

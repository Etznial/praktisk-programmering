using System;
using static System.Math;
using static System.Double;

public static partial class quad{

public static double adapt3o(
	Func<double,double> f,double a,double b,
	double acc=1e-3, double eps=1e-3, int limit=99,
	double f2=NaN)
{
	double f1=f(a+(b-a)/6), f3=f(a+5*(b-a)/6);
	if( IsNaN(f2) ) f2=f(a+3*(b-a)/6);
	double Q=(3*f1+2*f2+3*f3)/8*(b-a);
	double q=(f1+f2+f3)/3*(b-a);
	double err=Abs(Q-q)/2;
	if(limit==0){
		Console.Error.WriteLine($"adapt: limit reached: a={a} b={b}");
		return Q;
		}
	if(err<acc+eps*Abs(Q)){
		return Q;
		}
	else{
double Q1=adapt3o(f,a          ,a+(b-a)/3  ,acc/Sqrt(3),eps,limit-1,f1);
double Q2=adapt3o(f,a+(b-a)/3  ,a+2*(b-a)/3,acc/Sqrt(3),eps,limit-1,f2);
double Q3=adapt3o(f,a+2*(b-a)/3,b          ,acc/Sqrt(3),eps,limit-1,f3);
		return Q1+Q2+Q3;
		}
}
}

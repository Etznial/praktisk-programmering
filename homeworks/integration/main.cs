using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;

class main{
	static double integrate(
			Func<double,double> f, 
			double a, 
			double b,
			double delta=0.001, 
			double eps=0.001, 
			double f2=double.NaN, 
			double f3=double.NaN)
	{
		double h=b-a;
		if(double.IsNaN(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); } // first call, no points to reuse
		double f1=f(a+h/6), f4=f(a+5*h/6);
		double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); // lower order rule
		double err = Abs(Q-q);
		if (err <= delta+eps*Abs(Q)) return Q;
		else return integrate(f,  a     , (a+b)/2, delta/Sqrt(2), eps, f1, f2)+
			    integrate(f, (a+b)/2,    b   , delta/Sqrt(2), eps, f3, f4);
	}

 	static double erff1(double x){
		return Exp(-Pow(x,2));
	}
	static double erff2(double z, double t){
		return Exp(-Pow(z+(1-t)/t,2))/t/t;
	}

	static double erf(double z, double delta=0.001, double eps=0.001){
		if(z<0){ 
			return -erf(-z);
		}
		if(0<=z && z<=1){ 
			return 2.0/Sqrt(PI)*integrate(erff1,0,z,delta,eps);
		}
		else{
			Func<double,double> f3 = t => erff2(z,t);
			return 1.0-2.0/Sqrt(PI)*integrate(f3,0,1,delta,eps);

		}

	}
	
	
	static double f1(double x){
		return Sqrt(x);
	}
		
	static double f2(double x){
		return 1.0/Sqrt(x);
	}

	static double f3(double x){
		return 4.0*Sqrt(1-Pow(x,2));
	}

	static double f4(double x){
		return Log(x)/Sqrt(x);
	}


	public static void Main(){
		// A	
		WriteLine($"test integrate for sqrt(x)\t\tfrom 0 to 1, should give 2/3:\t{      approx(2.0/3,integrate(f1,0,1,1e-10,1e-10))}");
		WriteLine($"test integrate for 1/sqrt(x)\t\tfrom 0 to 1, should give 2:\t{      approx(2.0  ,integrate(f2,0,1,1e-10,1e-10))}");
		WriteLine($"test integrate for 4*sqrt(1-x^2)\tfrom 0 to 1, should give π:\t{  approx(PI   ,integrate(f3,0,1,1e-10,1e-10))}");
		WriteLine($"test integrate for Log(x)/sqrt(x)\tfrom 0 to 1, should give -4:\t{approx(-4   ,integrate(f4,0,1,1e-10,1e-10))}");
		
		WriteLine($"\nTest of approximation of erf from plots with table values and new and improved erf func with integrater with delta=1e-9 and eps=1e-9");
		WriteLine($"for z=1 \ntable value:\t0.842700793\napproximation:\t{sfuns.erf(1)}\nimproved erf:\t{erf(1,1e-9,1e-9)}\n");
		WriteLine($"for z=2 \ntable value:\t0.995322265\napproximation:\t{sfuns.erf(2)}\nimproved erf:\t{erf(2,1e-9,1e-9)}\n");
		WriteLine($"for z=3 \ntable value:\t0.999977910\napproximation:\t{sfuns.erf(3)}\nimproved erf:\t{erf(3,1e-9,1e-9)}\n");
		WriteLine("we can see from the above values that the impoved erf is more accurate than the approximate erf with a delta and eps of 1e-9");
		 
		
		
		string toWrite="";
		for(double i=0;i<=3.1;i+=0.1){
			toWrite+=$"{i}\t{erf(i)}\n";
		}
		File.WriteAllText("erf.data",toWrite);
		


	}	
}// class


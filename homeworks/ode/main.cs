using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;

class main{
	public static (vector, vector) rkstep12(
		Func<double,vector,vector> f, 	/* the f from dy/dx=f(x,y) */
		double x, 			/* the current value of the variable */
		vector y,			/* the current value y(x) of the sought function */
		double h			/* the step to be taken */
		)
	{
		vector k0 = f(x,y);		 /* embedded lower order formula (Euler) */
		vector k1 = f(x+h/2,y+k0*(h/2));  /* higher order formula (midpoint) */
		vector yh = y+k1*h;              /* y(x+h) estimate */
		vector er = (k1-k0)*h;           /* error estimate */
		return (yh,er);
	}
	public static (genlist<double>,genlist<vector>) driver(
			Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
			double a,			/* the start-point a */	
			vector ya, 			/* y(a) */
			double b,			/* the end-point of the integration */
			double h=0.01,			/* initial step-size */
			double acc=0.01,		/* absolute accuracy goal */
			double eps=0.01			/* relative accuracy goal */
		){
		if(a>b) throw new ArgumentException("driver: a>b");
		double x=a; vector y=ya.copy();
		var xlist=new genlist<double>(); xlist.add(x);
		var ylist=new genlist<vector>(); ylist.add(y);
		do{
			if(x>=b) return (xlist,ylist); /* job done */
			if(x+h>b) h=b-x;
			var (yh,erv) = rkstep12(f,x,y,h);
			double tol = (acc+eps*yh.norm()) * Sqrt(h/(b-a));
			double err = erv.norm();
			if(err<=tol){ // accept step
				x+=h; y=yh;
				xlist.add(x);
				ylist.add(y);
			}
			h *= Min( Pow(tol/err,0.25)*0.95 , 2); // reajust stepsize
			} while(true);
		} // driver


	public static vector f(double x, vector y){
		var y_dif=new vector(2);
		y_dif[0]= y[1];
		y_dif[1]=-y[0];
		return y_dif;
	}
	
	public static vector theta2(double t, vector y){
		double b=0.25, c=5;
		var y_diff=new vector(2);
		y_diff[0]=y[1]; 
		y_diff[1]=-b*y[1]-c*Sin(y[0]);
		return y_diff;
	}

	//exam
	//Finite difference formula for Hessian matrix
	public static matrix 
	

	public static void Main(){
		// rk test
		string toWriterk="";
		double a=0;
		double b=100;
		double h=0.5;
		vector y=new vector(2); y[0]=0; y[1]=1;
		for(int i = 0; i<b; i++){
			(vector yh,_) = rkstep12(f,h*i,y,h);
			toWriterk+=$"{i*h}\t{y[0]}\t{y[1]}\n";
			y=yh;
		}
		File.WriteAllText("rk.data",toWriterk);
		// Driver test
		vector ya=new vector(2); ya[0]=0; ya[1]=1;
		(genlist<double> xlist,genlist<vector> ylist)=driver(f,a,ya,b);
		string toWrite="";
		
		for(int i=0; i<xlist.size; i++){
			toWrite+=$"{xlist[i]}\t{ylist[i][0]}\n";
		}
		File.WriteAllText("driver.data",toWrite);
		//examples from scipy
		a=0;
		b=10;
		vector yas=new vector(2); yas[0]=PI-0.1; yas[1]=0.0;
		(genlist<double> xlist2,genlist<vector> ylist2)=driver(theta2,a,yas,b);
		toWrite="";
		string toWrite1="";
		for(int i=0; i<xlist2.size; i++){
			toWrite+=$"{xlist2[i]}\t{ylist2[i][0]}\n";
			toWrite1+=$"{xlist2[i]}\t{ylist2[i][1]}\n";

		}
		File.WriteAllText("pen.data",toWrite);
		File.WriteAllText("pen1.data",toWrite1);
		
	}	
}// class


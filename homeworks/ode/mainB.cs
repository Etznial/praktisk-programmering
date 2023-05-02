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
	public static vector driver(
			Func<double,vector,vector> f,	/* the f from dy/dx=f(x,y) */
			double a,			/* the start-point a */	
			vector ya, 			/* y(a) */
			double b,			/* the end-point of the integration */
			double h=0.01,			/* initial step-size */
			double acc=1e-2,		/* absolute accuracy goal */
			double eps=1e-2,			/* relative accuracy goal */
			genlist<double> xlist=null, genlist<vector> ylist=null
			){
		if(a>b) throw new ArgumentException("driver: a>b");
		if(xlist==null && ylist!=null) throw new ArgumentException("driver xlist=null while ylist!=null"); 
		if(xlist!=null && ylist==null) throw new ArgumentException("driver xlist!=null while ylist=null"); 
		double x=a; vector y=ya.copy();
		int xi=0;
		if(xlist[0]==a){
			ylist.add(y);
			xi++;
		}
		bool onPoint=false;
		do{
			if(x>=b){
				if(xlist[xi]==b) ylist.add(y);
				return y; /* job done */
			} 
			if(x+h>b) h=b-x;
			if(xi<xlist.size && xlist[xi]<x+h){
				h=xlist[xi]-x;
				onPoint=true;
			}
			var (yh,erv) = rkstep12(f,x,y,h);
			vector tol = new vector(y.size);
			for(int i=0;i<y.size;i++) tol[i]=(acc+eps*Abs(yh[i]))*Sqrt(h/(b-a));
			bool ok=true;
			for(int i=0;i<y.size;i++) if(!(erv[i]<tol[i])) ok=false;
			if(ok){ // accept step
				if(onPoint==true){
					ylist.add(y);
					onPoint=false;
					xi++;
				}
			       	x+=h; y=yh;

			}
			double factor = tol[0]/Abs(erv[0]);
			for(int i=1;i<y.size;i++) factor=Min(factor,tol[i]/Abs(erv[i]));
			h *= Min( Pow(factor,0.25)*0.95 ,2);
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

	public static vector u(double t, vector y){
		double eps = 0;
		var y_diff=new vector(2);
		y_diff[0]=y[1]; 
		y_diff[1]=1+eps*Pow(y[0],2)-y[0];
		return y_diff;
	}
	public static vector u1(double t, vector y){
		double eps = 0.01;
		var y_diff=new vector(2);
		y_diff[0]=y[1]; 
		y_diff[1]=1+eps*Pow(y[0],2)-y[0];
		return y_diff;
	}

	public static void Main(){
		double a=0;
		double b=10;
		genlist<vector> ylist = new genlist<vector>();
		genlist<double> xlist = new genlist<double>();
		for(int i=0;i<11;i++){
			xlist.add(i);
		}

		vector ya=new vector(2); ya[0]=PI-0.1; ya[1]=0.0;
		vector y=driver(theta2,a,ya,b, xlist:xlist, ylist:ylist);
		string toWrite="";
		for(int i=0; i<xlist.size; i++){
			toWrite+=$"{xlist[i]}\t{ylist[i][0]}\n";

		}
		File.WriteAllText("penB.data",toWrite);
		
		double a_phi = 0;
		double b_phi = 2*PI*10;
		genlist<vector> ylist1 = new genlist<vector>();
		genlist<vector> ylist2 = new genlist<vector>();
		genlist<vector> ylist3 = new genlist<vector>();
		genlist<double> xlist1 = new genlist<double>();
		double res = 2000;
		for(int i=0;i<=res;i++){
			xlist1.add(i/res*b_phi);
		}
		

		vector ya1 = new vector(2); ya1[0]=1; ya1[1]=0;
		vector y1=driver(u,a_phi,ya1,b_phi, xlist:xlist1, ylist:ylist1);
		vector ya2 = new vector(2); ya2[0]=1; ya2[1]=-0.5;
		vector y2=driver(u,a_phi,ya2,b_phi, xlist:xlist1, ylist:ylist2);
		vector ya3 = new vector(2); ya3[0]=1; ya3[1]=-0.5; // use with eps = 0.01
		vector y3=driver(u1,a_phi,ya3,b_phi, xlist:xlist1, ylist:ylist3);
		string toWrite1="";
		string toWrite2="";
		string toWrite3="";
		for(int i=0;i<xlist1.size;i++){
			toWrite1+=$"{xlist1[i]}\t{ylist1[i][0]}\n";
			toWrite2+=$"{xlist1[i]}\t{ylist2[i][0]}\n";
			toWrite3+=$"{xlist1[i]}\t{ylist3[i][0]}\n";
		}
		File.WriteAllText("u1.data",toWrite1);
		File.WriteAllText("u2.data",toWrite2);
		File.WriteAllText("u3.data",toWrite3);
	}	
}// class


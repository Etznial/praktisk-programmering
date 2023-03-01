using System;
using static System.Math;
using static System.Console;
static class main{
static double pow(this double x,int n){return Pow(x,n);}
static double pow(this double x,double n){return Pow(x,n);}
public static void Main(){
Func<vector,double> f;

f = t => {
	double A=10,sum=A*t.size;
	for(int i=0;i<t.size;i++) sum+=t[i]*t[i]-A*Cos(2*PI*t[i]);
	return sum;
	};

f = t => {
	double A=100,sum=0;
	for(int i=0;i<t.size-1;i++)
		sum+=A*(t[i+1]-t[i]*t[i]).pow(2)+(1-t[i])*(1-t[i]);
	return sum;
	};

var s=new vector("0 0");
vector r=swarm.run(f,start:s,step:5.12,seconds:1);
r.print("r=");
WriteLine($"f(r)={f(r)}");
}
}//class

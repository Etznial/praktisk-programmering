using System;
using System.Diagnostics;
using static System.Console;

public class halton{
int[] bs = new int[]{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
int n=0,s=0,d=1;
	public static double corput(int n, int b){
		double q=0,bk=1.0/b;
		while(n>0){q+=(n%b)*bk;n/=b;bk/=b;}
		return q;
	}
	public halton(int dim=1,int shift=0){d=dim;s=shift;}
	public vector next(){
		n++;
		vector x=new vector(d);
		for(int i=0;i<d;i++)x[i]=corput(n,bs[i+s]);
		return x;
	}
}

public static class swarm{
public static System.Random rnd=new Random(1);

public static vector randomvec(vector a,vector b){
	var x=new vector(a.size);
	for(int i=0;i<a.size;i++)x[i]=a[i]+(b[i]-a[i])*rnd.NextDouble();
	return x;
	}

public static vector run
(Func<vector,double>f,vector start, double step,int seconds=1,double w=0.723)
{
vector a=start.copy();
vector b=start.copy();
for(int i=0;i<start.size;i++){
	a[i]-=step;
	b[i]+=step;
	}
rnd=new Random(2);
int dim=a.size;
int N=8*dim;
vector[] x=new vector[N], v=new vector[N], p=new vector[N];
var fp=new double[N];
var g=randomvec(a,b);
double fg=f(g);
for(int i=0;i<N;i++){
	x[i]=randomvec(a,b);
	p[i]=x[i].copy();
	fp[i]=f(p[i]);
	if(fp[i]<fg){g=p[i].copy();fg=fp[i];}
	v[i]=randomvec(a-b,b-a)/2;
	}
var time=DateTime.Now;
do{
	for(int i=0;i<N;i++){
		v[i]=w*v[i]
			+rnd.NextDouble()*(p[i]-x[i])
			+rnd.NextDouble()*(g-x[i]);
		x[i]+=v[i];
		var fxi=f(x[i]);
		if(fxi<fp[i]){fp[i]=fxi;p[i]=x[i].copy();}
		if(fxi<fg)   {fg=fxi;   g=x[i].copy();}
	}
}while((DateTime.Now-time).TotalSeconds < seconds);
return g;	
}

}//class

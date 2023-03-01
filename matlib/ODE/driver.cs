using System;
using static System.Console;
using static System.Math;
using System.Collections.Generic;
public partial class ODE{

public static
Func<Func<double,vector,vector>,double,vector,double,(vector,vector)>
stepper = rkstep23;
//stepper = rkstep12;

public static vector driver(
	Func<double,vector,vector> F, /* equation */
	double a, vector ya, double b, 
	double acc=0.01, double eps=0.01, double h=0, 
	List<double> xlist=null, List<vector> ylist=null,
	int limit=999
	)
{// ODE driver
int nsteps=0;
if(h==0)h=(b-a)/64;
if(xlist!=null) {xlist.Clear(); xlist.Add(a);}
if(ylist!=null) {ylist.Clear(); ylist.Add(ya);}
double x=a;
vector y=ya, tol=new vector(y.size);
do{
	if(x>=b) return y;
	if(nsteps>limit) throw new Exception($"driver: nsteps>{limit}\n");
	if(x+h>b) h=b-x;
	var (yh,erv) = stepper(F,x,y,h);

	for(int i=0;i<tol.size;i++) tol[i]=Max(acc,Abs(yh[i])*eps)*Sqrt(h/(b-a));
	bool ok=true; for(int i=0;i<tol.size;i++) if(Abs(erv[i])>tol[i])ok=false;

//	double tol = Max(acc,yh.norm()*eps) * Sqrt(h/(b-a));
//	double err = erv.norm();
//	bool ok = (err<=tol);

	if(ok){
		x+=h; y=yh; nsteps++;
		if(xlist!=null) xlist.Add(x);
		if(ylist!=null) ylist.Add(y);
		}

//	if(err==0) h *= 2;
//	else h *= Min( Pow(tol/err,0.25)*0.95 , 2);

	double factor = tol[0]/Abs(erv[0]);
	for(int i=1;i<tol.size;i++) factor=Min(factor,tol[i]/Abs(erv[i]));
	h *= Min( Pow(factor,0.25)*0.95 ,2);

	}while(true);

}//driver

}//class

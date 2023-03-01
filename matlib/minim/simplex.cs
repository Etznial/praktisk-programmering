using System;
using static System.Console;
public static class simplex{

public static double size(vector[] p){
	double s=0;
	for(int i=1;i<p.Length;i++){ s=Math.Max(s,(p[0]-p[i]).norm()); }
	return s;
	}

public static vector downhill
(System.Func<vector,double>F,vector x,double step=1.0/64,double minsize=1e-3){

vector[] p=new vector[x.size+1];
double[] f=new double[x.size+1];
p[x.size]=x.copy();
f[x.size]=F(p[x.size]);
for(int i=0;i<x.size;i++){
	x[i]+=step;
	p[i]=x.copy();
	f[i]=F(p[i]);
	x[i]-=step;
	}
int hi=0,lo=0,nsteps=0;
vector pnew; double fnew;

while(size(p)>minsize && ++nsteps<999){
	hi=0;lo=0;
	double fhi=f[hi],flo=f[lo];
	for(int k=1;k<p.Length;k++){
		if(f[k]>fhi){fhi=f[k];hi=k;}
		if(f[k]<flo){flo=f[k];lo=k;}
		}
	vector pce=new vector(p[0].size);
	for(int i=0;i<p.Length;i++) if(i!=hi)pce+=p[i];
	pce/=pce.size;
	pnew=3*pce-2*p[hi];
// expansion
	fnew=F(pnew);
	if(fnew<flo){
		p[hi]=pnew;
		f[hi]=fnew;
		continue;
		}
// reflection
	pnew=2*pce-p[hi];
	fnew=F(pnew);
	if(fnew<fhi){
		p[hi]=pnew;
		f[hi]=fnew;
		continue;
		}
// contraction
	pnew=(pce+p[hi])/2;
	fnew=F(pnew);
	if(fnew<fhi){
		p[hi]=pnew;
		f[hi]=fnew;
		continue;
		}
// reduction
	for(int i=0;i<p.Length;i++)
		if(i!=lo){
			p[i]=(p[i]+p[lo])/2;
			f[i]=F(p[i]);
		}
	}//while
return p[lo];
}//downhill
}//class simplex

using System;
using static System.Math;
public static partial class roots{

public static readonly double DELTA = Pow(2,-26);

public static matrix jacobian
(Func<vector,vector> f,vector x,vector fx=null,vector dx=null){
	int n=x.size;
//	if(dx==null) dx = x.map(t => Max(Abs(t),1)*DELTA);
	if(dx==null) dx = x.map(t => Abs(t)*DELTA);
	matrix J=new matrix(n,n);
	if(fx == null)fx=f(x);
	for(int j=0;j<n;j++){
		x[j]+=dx[j];
		vector df=f(x)-fx;
		for(int i=0;i<n;i++) J[i,j]=df[i]/dx[j];
		x[j]-=dx[j];
		}
	return J;
}

}//roots

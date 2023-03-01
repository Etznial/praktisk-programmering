using System;
using static System.Math;
using static complex;
using static cmath;

public static partial class matlib{

public static void dfts
(int sign, int N, complex[] x, int ix, int stride, complex[] c, int ic){
	for(int k=0;k<N;k++){
		c[ic+k]=0;
		for(int n=0;n<N;n++)
			c[ic+k]+=x[ix+n*stride]*exp(sign*2*PI*I*n*k/N);
	}
}

public static void ffts
(int sign, int N, complex[] x, int ix, int stride, complex[] c, int ic){
	if(N==1) c[ic+0]=x[ix+0];
	else if(N%2==0){
		ffts(sign,N/2,x,ix+0     ,2*stride,c,ic+0  );
		ffts(sign,N/2,x,ix+stride,2*stride,c,ic+N/2);
		for(int k=0;k<N/2;k++){
			complex p=c[ic+k], q=exp(sign*2*PI*I*k/N)*c[ic+k+N/2];
			c[ic+k    ]=p+q;
			c[ic+k+N/2]=p-q;
			}
		}
	else dfts(sign,N,x,ix,stride,c,ic);
	}

public static complex[] fft(complex[] x){
	int N=x.Length;
	var c=new complex[N];
	ffts(-1,N,x,0,1,c,0);
	return c;
	}

public static complex[] ift(complex[] c){
	int N=c.Length;
	var x=new complex[N];
	ffts(+1,N,c,0,1,x,0);
	for(int i=0;i<N;i++)x[i]/=N;
	return x;
	}

}//class

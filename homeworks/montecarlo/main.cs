using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;

class main{
	public static (double,double) plainmc(Func<vector,double> f,vector a,vector b,int N){
		int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
		double sum=0,sum2=0;
		var x=new vector(dim);
		var rnd=new Random();
		for(int i=0;i<N;i++){
			for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
			double fx=f(x); sum+=fx; sum2+=fx*fx;
		}
		double mean=sum/N, sigma=Sqrt(sum2/N-mean*mean);
		var result=(mean*V,sigma*V/Sqrt(N));
		return result;
	}

	public static (double,double) qrmc(Func<vector,double> f, vector a, vector b, int N){
		int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
		double sum0=0,sum1=0;
		var x0=new vector(dim);
		var x1=new vector(dim);
		var h0=new vector(dim);
		var h1=new vector(dim);
		for(int i=0;i<N;i++){
			halton(i,dim,h0,0);
			for(int k=0;k<dim;k++)x0[k]=a[k]+h0[k]*(b[k]-a[k]);
			halton(i,dim,h1,dim);
			for(int k=0;k<dim;k++)x1[k]=a[k]+h1[k]*(b[k]-a[k]);
			double fx0=f(x0); sum0+=fx0; 
			double fx1=f(x1); sum1+=fx1;
		}
		double mean0=sum0/N;
		double mean1=sum1/N;
		double inte0=mean0*V;
		double inte1=mean1*V;
		double tempSum=0;
		double mean01=(inte0+inte1)/2;
		tempSum+=Pow(inte0-mean01,2);
		tempSum+=Pow(inte1-mean01,2);
		double sigma=Sqrt(1.0/2*tempSum);
		var result=(mean01,sigma/Sqrt(2));
		return result;
	}

	public static double corput(int n, int b){
		double q=0;
		double bk=1.0/b;
		while(n>0){
			q+=(n%b)*bk; 
			n/=b; 
			bk/=b;
		}
		return q;
	}

	public static void halton(int n, int d, vector x, int k=0){
		int[] bs={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
		int maxd=bs.Length/d;
		Debug.Assert(d<=maxd,"dimension is too big, need larger basis"); 
		for(int i=0; i<d; i++) x[i]=corput(n,bs[i+k]);
	}
	
	public static void Main(){
		Func<vector, double> cossin = (v) => Cos(v[0])+Sin(v[1]);
		Func<vector, double> coscoscos = (v) => Pow(1-Cos(v[0])*Cos(v[1])*Cos(v[2]),-1)/Pow(PI,3);

		vector vbcs = new vector(-PI,-PI); // limits for integral of cossin
		vector vecs = new vector( PI, PI); // limits for integral of cossin
		vector vbccc = new vector(0,0,0); // limits for integral of coscoscos
		vector veccc = new vector(PI,PI,PI); // limits for integral of coscoscos
		int N = (int)1e6;


		var testcossin = plainmc(cossin,vbcs,vecs,N);
		WriteLine("====================[A]====================");
		WriteLine("Test integral of cos x + sin y with limits -pi to pi for both x and y, the integral is calculated to 0 analytically");
		WriteLine($"result of plainmc: {testcossin}");
		WriteLine("Test integral of [1-cos(x)cos(y)cos(z)]**(-1)/PI**3, result should approximatly be:");
		WriteLine("1.3932039296856768591842462603255");
		WriteLine($"result of plainmc: {plainmc(coscoscos,vbccc,veccc,N)}");
		WriteLine("for referance see plot Testplainmc.svg");
		//Time to make graph
		int beginIt = (int)50;
		int endIt = (int)1e5;
		string toWrite="";
		for(int n=beginIt; n<endIt; n+=beginIt+n/10){
			(double integral, double error) = plainmc(cossin,vbcs,vecs,n);
			toWrite+=$"{n}\t{error}\n";
		}
		File.WriteAllText("coscosIt.data",toWrite);
		
		WriteLine("");
		WriteLine("====================[B]====================");
		WriteLine("Estimate the error by using two different sequences");
		//Time to make graph
		toWrite="";
		for(int n=beginIt; n<endIt; n+=beginIt+n/10){
			(double integral, double error) = qrmc(cossin,vbcs,vecs,n);
			toWrite+=$"{n}\t{error}\n";
		}
		File.WriteAllText("QRMCcoscosIt.data",toWrite);
		WriteLine("The compareson of the scaling of the errors of plainmc and qrmc is done in the plot Testqrmc.svg");



	}	
}// class


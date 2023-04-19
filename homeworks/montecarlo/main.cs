using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

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

	//make 2D function
	
	
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
		
		//Time to make graph
		int beginIt = (int)1e1;
		int endIt = (int)1e4;
		string toWrite="";
		for(int n=beginIt; n<endIt; n+=10){
			(double integral, double error) = plainmc(cossin,vbcs,vecs,n);
			toWrite+=$"{n}\t{error}\n";
		}

		File.WriteAllText("coscosIt.data",toWrite);
		



	}	
}// class


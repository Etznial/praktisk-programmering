using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class main{
	
	static double polinterp(int n, double[] x, double[] y, double z){
		double s=0,p;
		for(int i=0;i<n;i++) {
			p=1;
			for(int k=0;k<n;k++){
				if(k!=i){
					p*=(z-x[k])/(x[i]-x[k]);
				}
			}
			s+=y[i]*p;
		}
		return s;
	}
	public static int binsearch(double[] x, double z){/* locates the interval for z by bisection */ 
				if(!(x[0]<=z && z<=x[x.Length-1])) throw new Exception("binsearch: bad z");
				int i=0, j=x.Length-1;
				while(j-i>1){
				int mid=(i+j)/2;
				if(z>x[mid]) i=mid; else j=mid;
		}
		return i;
	}

	public static double linterp(double[] x, double[] y, double z){
		        int i=binsearch(x,z);
			double dx=x[i+1]-x[i]; 
			if(!(dx>0)) throw new Exception("dx needs to be dx>0");
			double dy=y[i+1]-y[i];
			return y[i]+dy/dx*(z-x[i]);
	}

	public static double linterpInteg(double[] x, double[] y, double z){
		double sum=0;
		int i = binsearch(x,z);
		
		for(int j=0; j<i; j++){
			double a=(y[j+1]-y[j])/(x[j+1]-x[j]);
			double b=y[j];
			sum+=0.5*a*Pow((x[j+1]-x[j]),2)+b*(x[j+1]-x[j]);
		}
		double ai=(y[i+1]-y[i])/(x[i+1]-x[i]);
		double bi=y[i];
		sum+=0.5*ai*Pow((z-x[i]),2)+bi*(z-x[i]);
		return sum;
	}

	class qspline {
		vector x,y,p,c;
		public qspline(vector xs,vector ys){
		x=xs.copy(); y=ys.copy(); 
		int n=x.size;
		p = new vector(n-1);
		c = new vector(n-1);
		for(int i=0; i<n-1; i++){
			p[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
		}
		c[0]=0;
		for(int i=0; i<n-2; i++){
			double dx =x[i+1]-x[i];
			double dx1=x[i+2]-x[i+1];
			c[i+1]=1/dx1*(p[i+1]-p[i]-c[i]*dx);
		}
		
		c[n-2]=c[n-2]*0.5;
		for(int i=n-3; i>=0; i--){
			double dx =x[i+1]-x[i];
			double dx1=x[i+2]-x[i+1];
			c[i]=1/dx*(p[i+1]-p[i]-c[i+1]*dx1);
		}
		
		
		}
		public double evaluate(double z){/* evaluate the spline */
			int i=binsearch(x,z);
			double s=y[i]+p[i]*(z-x[i])+c[i]*(z-x[i])*(z-x[i+1]);
			return s;
		}
		public double derivative(double z){/* evaluate the derivative */
			int i=binsearch(x,z);
			double bi=p[i]-c[i]*(x[i+1]-x[i]);
			return bi+2*c[i]*(z-x[i]);
		}
		public double integral(double z){/* evaluate the integral */
			int i=binsearch(x,z);
			double sum=0;
			for(int j=0; j<i; j++){
				double bj=p[j]-c[j]*(x[j+1]-x[j]);
				sum+=y[j]*(x[j+1]-x[j])+0.5*bj*Pow(x[j+1]-x[j],2)+1.0/3*c[j]*Pow(x[j+1]-x[j],3);
			}
			double bi=p[i]-c[i]*(x[i+1]-x[i]);
			sum+=y[i]*(z-x[i])+0.5*bi*Pow(z-x[i],2)+1.0/3*c[i]*Pow(z-x[i],3);
			return sum;
		}
		public double getbi(int i){	
			double bi=p[i]-c[i]*(x[i+1]-x[i]);
			return bi;
		}
		public double getci(int i){
			return c[i];
		}
	}


	public static void Main(){
		//make data
		Random rand = new Random(2);
		//double maxValue= 10;
		//double minValue=-10;
		int n = 20;
		string toWrite="";
		double[] xs = new double[n];
		double[] ys = new double[n];
		double yVal = rand.NextDouble();
		for(int i=0; i<n; i++){
			yVal += rand.NextDouble()-0.5;
			toWrite+=$"{i}\t{yVal}\n";
			xs[i] = i;
			ys[i] = yVal;
		}	
		File.WriteAllText("test.data",toWrite);
		File.WriteAllText("testQS.data",toWrite);
		//Fit
		toWrite="";
		double res = 1000;
		for(int i=0; i<res+1; i++){
			double z = i*((xs[n-1]-xs[0])/res);
			toWrite+=$"{z}\t{linterp(xs,ys,z)}\n";
		}
		File.WriteAllText("linSpline.data",toWrite);
		// integral
		double zEnd = xs[xs.Length-1];
		WriteLine($"Integralsum of lin spline: {linterpInteg(xs,ys,zEnd)}");
		double zRun = 0;
		double step = (xs[xs.Length-1]-xs[0])/res;
		toWrite="";
		while(zEnd>zRun){
			toWrite+=$"{zRun}\t{linterpInteg(xs,ys,zRun)}\n";
			zRun+=step;
		}
		File.WriteAllText("linIntegral.data",toWrite);
		// B
		//Fit
		toWrite="";
		res = 1000;
		var qspline = new qspline(xs,ys);
		for(int i=0; i<res+1; i++){
			double z = i*((xs[n-1]-xs[0])/res);
			toWrite+=$"{z}\t{qspline.evaluate(z)}\n";
		}
		File.WriteAllText("qspline.data",toWrite);
		
		res =1000;
		zEnd = xs[xs.Length-1];
		zRun = 0;
		step = (xs[xs.Length-1]-xs[0])/res;
		toWrite="";
		string toWriteD="";
		while(zEnd>zRun){
			toWrite+=$"{zRun}\t{qspline.integral(zRun)}\n";
			toWriteD+=$"{zRun}\t{qspline.derivative(zRun)}\n";
			zRun+=step;
		}
		File.WriteAllText("qsplineIntegral.data",toWrite);
		File.WriteAllText("qsplineDerivative.data",toWriteD);
			
		WriteLine("==========[B]==========");
		WriteLine("Calcullate manually the parameters {bi, ci} of the corresponding quadratic-splines, and compare the results with your quadratic-spline program.");
		vector x = new vector(5);
		for(int i=0; i<x.size; i++){
			x[i]=i+1;
		}
		vector y1 = new vector(5);
		vector y2 = new vector(5);
		vector y3 = new vector(5);
		for(int i=0; i<x.size; i++){
			y1[i]=1;
		}
		for(int i=0; i<x.size; i++){
			y2[i]=i+1;
		}
		for(int i=0; i<x.size; i++){
			y3[i]=Pow(i+1,2);
		}
		var qspliney1 = new qspline(x,y1);
		var qspliney2 = new qspline(x,y2);
		var qspliney3 = new qspline(x,y3);
		
		vector bi1 = new vector(4);
		vector bi2 = new vector(4);
		vector bi3 = new vector(4);
		
		for(int i=0; i<bi1.size; i++){
			bi1[i]=qspliney1.getbi(i);
		}
		for(int i=0; i<bi2.size; i++){
			bi2[i]=qspliney2.getbi(i);
		}
		for(int i=0; i<bi3.size; i++){
			bi3[i]=qspliney3.getbi(i);
		}
		
		vector ci1 = new vector(4);
		vector ci2 = new vector(4);
		vector ci3 = new vector(4);

		for(int i=0; i<ci1.size; i++){
			ci1[i]=qspliney1.getci(i);
		}
		for(int i=0; i<ci2.size; i++){
			ci2[i]=qspliney2.getci(i);
		}
		for(int i=0; i<ci3.size; i++){
			ci3[i]=qspliney3.getci(i);
		}
		
		WriteLine("x:");
		x.print();
		WriteLine("test with y1");
		WriteLine("y1:");
		y1.print();
		WriteLine("Calculated manually");
		WriteLine("bi1");
		WriteLine("\t0\t0\t0\t0");
		WriteLine("ci1");
		WriteLine("\t0\t0\t0\t0");
		WriteLine("");
		WriteLine("Calculated with program");
		WriteLine("bi1");
		bi1.print();
		WriteLine("ci1");
		ci1.print();
		WriteLine("");
		WriteLine("test with y2");
		WriteLine("y2:");
		y2.print();
		WriteLine("Calculated manually");
		WriteLine("bi2");
		WriteLine("\t1\t1\t1\t1");
		WriteLine("ci2");
		WriteLine("\t0\t0\t0\t0");
		WriteLine("");
		WriteLine("Calculated with program");
		WriteLine("bi2");
		bi2.print();
		WriteLine("ci2");
		ci2.print();
		WriteLine("");
		WriteLine("test with y3");
		WriteLine("y3:");
		y3.print();
		WriteLine("Calculated manually");
		WriteLine("bi3");
		WriteLine("\t2\t4\t6\t8");
		WriteLine("ci3");
		WriteLine("\t1\t1\t1\t1");
		WriteLine("");
		WriteLine("Calculated with program");
		WriteLine("bi3");
		bi3.print();
		WriteLine("ci3");
		ci3.print();
	
	
	}
}// class

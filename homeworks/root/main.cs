using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;

class main{
	public static int it=0;
	static vector newton(Func<vector,vector> f, vector x, double eps=1e-2){
		vector cx=x.copy();
		vector fcx = f(cx);
		cx.print("cx=");
		fcx.print("fcx=");
		int m = fcx.size; // rows
		int n = cx.size; // columbs
		if(n!=m) throw new ArgumentException($"n must be equal to m to have a square matrix, they are atm n={n} and m={m}");
		matrix J = new matrix(n,m);
		vector step = new vector(n);
		while(fcx.norm()>eps){
			for(int i=0;i<m;i++){ // it rows
				for(int k=0;k<n;k++){ // it columbs
					//if(Abs(x[i]) < Pow(2,-26)) x[i] = Pow(2,-24);
					double dxk = Abs(cx[k])*Pow(2,-26);
					cx[k]+=dxk;
					vector ff=f(cx);
					J[i,k]=(ff[i]-fcx[i])/dxk;
					cx[k]-=dxk;
				}	
			}
			step = QRGS.solveG(J,-fcx);
			step.print("step=");
			var mismatch=J*step+fcx;
			mismatch.print("skal være nul");
			double lambda = 1;
			vector z=cx+step*lambda;	
			vector fz=f(z);
			do{
				if(fz.norm()<(1-lambda/2)*fcx.norm())break;
				else{
					WriteLine("backtracking...");
					lambda/=2;
					z=cx+step*lambda;	
					fz=f(z);
					//WriteLine($"fz.norm(): {fz.norm()}");
					//WriteLine($"lambda: {lambda}");
					//WriteLine($"(1-lambda/2)*fcx.norm(): {(1-lambda/2)*fcx.norm()}");
				}
			}while( lambda>Pow(2,-23) );
			
			fcx=fz;
			cx=z;
			z.print("z=");
			fz.print($"|fz|={fz.norm()} fz=");
		}
		return cx;
	}
	/*
	static vector newtonTest(Func<vector,vector> f, vector x, double eps=1e-2){
		int n = x.size;
		int m = f(x).size;
		matrix J = new matrix(n,m);
		vector step = new vector(n);
		double dxk;
		double lambda;
		while(f(x).norm()>eps){
			//WriteLine($"{f(x).norm()}\t{eps}");
			// make J
			for(int i=0;i<m;i++){
				for(int k=0;k<n;k++){
					dxk = Abs(x[k])*Pow(2,-26);
					vector xdx = x.copy(); xdx[k]+=dxk;
					J[i,k]=(f(xdx)[i]-f(x)[i])/dxk;
				}
			}
			step = QRGS.solveG(J,-f(x));
			lambda = 1;
			while(f(x+step).norm() > (1.0-lambda*0.5)*f(x).norm() && lambda > 1.0/32.0) lambda/=2;
			x += lambda*step;
		}
		return x;
	}
	*/
	
	static public vector f(vector v){
		it++;
		var res = new vector(2);
		res[0]=-2*(1-v[0])-4*(v[1]-Pow(v[0],2))*v[0];
		res[1]=2*(v[1]-Pow(v[0],2));
		return res;
	}

	static public vector test1(vector v){
		it++;
		var res = new vector(v.size);
		for(int i=0;i<v.size;i++)res[i]=v[i]*v[i];
		return res;
	}
		

	public static void Main(){
		
		WriteLine("newton test");
		vector start;
		start = new vector(-5,7);
		it=0;
		newton(f,start).print();
		WriteLine($"number of times func has been called: {it}");
		start = new vector(1,1,1);
		it=0;
		newton(test1,start,eps:1e-5).print();
		WriteLine($"number of times func has been called: {it}");
		
		/*
		string toWrite1="";
		toWrite1+=$"{xlist2[i]}\t{ylist2[i][1]}\n";
		File.WriteAllText("pen1.data",toWrite1);
		*/
		/*
		WriteLine("newton test");
		vector start;
		start = new vector(-5,7);
		it=0;
		newtonTest(f,start).print();
		WriteLine($"number of times func has been called: {it}");
		start = new vector(1,1,1);
		it=0;
		newtonTest(test1,start,eps:1e-5).print();
		WriteLine($"number of times func has been called: {it}");
		*/
	}	
}// class


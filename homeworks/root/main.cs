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
		//cx.print("cx=");
		//fcx.print("fcx=");
		int m = fcx.size; // rows
		int n = cx.size; // columbs
		double lambda;
		double dxk;
		if(n!=m) throw new ArgumentException($"n must be equal to m to have a square matrix, they are atm n={n} and m={m}");
		matrix J = new matrix(n,m);
		vector step = new vector(n);
		while(fcx.norm()>eps){
			
			for(int i=0;i<m;i++){ // it rows
				for(int k=0;k<n;k++){ // it columbs
					//if(Abs(x[i]) < Pow(2,-26)) x[i] = Pow(2,-24);
					dxk = Abs(cx[k])*Pow(2,-26);
					cx[k]+=dxk;
					vector ff=f(cx);
					J[i,k]=(ff[i]-fcx[i])/dxk;
					cx[k]-=dxk;
				}	
			}
			
			step = QRGS.solveG(J,-fcx);
			//step.print("step=");
			//var mismatch=J*step+fcx;
			//mismatch.print("skal være nul");
			//J.print("Jacobian matrix: ");
			lambda = 1;
			vector z=cx+step*lambda;	
			vector fz=f(z);
			
			do{
				if(fz.norm()<(1-lambda/2)*fcx.norm())break;
				else{
					//WriteLine("backtracking...");
					lambda/=2;
					z=cx+step*lambda;	
					fz=f(z);
				}
			}while( lambda>Pow(2,-23) );
			fcx=fz;
			cx=z;
			//z.print("z=");
			//fz.print($"|fz|={fz.norm()} fz=");
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
			step.print("step=");
			var mismatch=J*step+f(x);
			mismatch.print("skal være nul");
			lambda = 1;
			while(f(x+lambda*step).norm() > (1.0-lambda*0.5)*f(x).norm() && lambda > 1.0/32.0) {
				WriteLine("backtracking...");
				WriteLine($"f(x+step).norm()={f(x+lambda*step).norm()}");
				lambda/=2;
			} 
			x += lambda*step;
		}
		return x;
	}
	*/
	/*
	static vector newtonTest(Func<vector,vector> f, vector x, double eps=1e-2){
		int m = x.size;
		int n = f(x).size;
		double delta_x;
		vector x1;
		double lambda;
		if(n!=m) WriteLine("Function vector and variable vector must be same size");
		matrix J = new matrix(n,m);
		while(f(x).norm() > eps){
			for(int i = 0; i<n; i++){
				if(Abs(x[i]) < Pow(2,-26)) x[i] = Pow(2,-24);
					delta_x = Abs(x[i])*Pow(2,-26);
					for(int k = 0; k<m; k++){
					  	x1 = x.copy(); x1[k] += delta_x;
						J[i,k] = (f(x1)[i] - f(x)[i])/delta_x;
					}
				}
				J.print("Jacobian matrix: ");
				vector d_x = QRGS.solveG(J,-f(x));
				d_x.print("step=");
				var mismatch=J*d_x+f(x);
				mismatch.print("skal være nul");
				lambda = 1;
				while(f(x+lambda*d_x).norm() > (1-0.5*lambda)*f(x).norm() && lambda > 1.0/32.0){
				       	lambda /= 2;
					WriteLine("backtracking...");
				}
				x += lambda*d_x;
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
	static public vector f2(vector v){
		it++;
		var res = new vector(2);
		res[0]=Abs(v[0])+Abs(v[1]);
		res[1]=Abs(v[0])+Abs(v[1]);
		return res;
	}


	static public vector test1(vector v){
		it++;
		var res = new vector(v.size);
		for(int i=0;i<v.size;i++)res[i]=v[i]*v[i];
		return res;
	}
		
	static public vector test2(vector v){
		it++;
		var res = new vector(2);
		res[0]=2*(v[0]+1);
		res[1]=2*v[1];
		return res;
	}

	static public vector test3(vector v){
		it++;
		var res = new vector(2);
		res[0]=2*v[0];
		res[1]=2*v[1];
		return res;
	}

	static public vector test4(vector v){
		it++;
		var res = new vector(2);
		res[0]=v[1];
		res[1]=v[0];
		return res;
	}

	static public vector test5(vector v){
		it++;
		var res = new vector(2);
		res[0]=2*(v[0]-1);
		res[1]=2*(v[1]+1);
		return res;
	}
	

	public static void Main(){
		// A
		// testing with tests 1 through 4
		WriteLine("=====[A]=====");
		int a = 5;
		double acc=1e-10;
		vector start = new vector(a,a);
		start.print($"starting parameters for all tests:");
		WriteLine($"acc={acc}\n");
		
		WriteLine("x**2+y**2");
		newton(test1, start,acc).print("test1 awnser:");
		WriteLine($"should have extremum at (0,0). The test function was called {it} times\n"); it=0;

		start = new vector(a,a);
		WriteLine("(x+1)**2+y**2");
		newton(test2, start,acc).print("test2 awnser:");
		WriteLine($"should have extremum at (-1,0). The test function was called {it} times\n"); it=0;

		start = new vector(a,a);
		WriteLine("test of x**2+y**2");
		newton(test3, start,acc).print("test3 awnser:");
		WriteLine($"should have extremum at (0,0). The test function was called {it} times\n"); it=0;
		
		start = new vector(a,a);
		WriteLine("test of x+y");
		newton(test4, start,acc).print("awnser:");
		WriteLine($"should have extremum at (0,0). The test function was called {it} times\n"); it=0;

		start = new vector(a,a);
		WriteLine("test of (x-1)**2+(y+1)**2");
		newton(test5, start,acc).print("awnser:");
		WriteLine($"should have extremum at (1,-1). The test function was called {it} times\n"); it=0;

		it=0;
		WriteLine("The Rosenbrock's valley function");
		start = new vector(10,-10);
		start.print("startinf parameters for the Rosenbrock's valley function");
		newton(f,start).print("awnser: ");
		WriteLine($"should have extremum at (1,1). The Rosenbrock's valley function was called {it} times\n"); it=0;
		WriteLine("The Rosenbrock's valley function is made to be hard for minimizers, because of it's sharp vally, which explains the absurd amount of calls to the function\n");
		//B	
		WriteLine("=====[B]=====");

	
	
	}	
}// class


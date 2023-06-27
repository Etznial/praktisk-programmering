using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;
using System.Collections.Generic;

class main{
	public static int it=0;
	public const double MEPS=2e-10;
	public static int step=0;
	public static vector qnewton(
		Func<vector,double> f, 	/* objective function */
		vector start, 		/* starting point */
		double acc = 1e-3	/* accuracy goal, on exit |gradient| should be < acc */){
		
		int n = start.size;
		
		double lambda;
		double uy;
		vector dx;
		vector s; /* s=lambda*dx */
		vector u;
		vector y;
		vector grad=gradient(f,start);
		vector gradOld;
		matrix dB;
		matrix B = new matrix(n,n);
		B.set_identity();
		for(int i=0; i<n; i++){
			if(Abs(start[i])<Pow(2,-26)) start[i]=Pow(2,-26);
		}
		while(grad.norm()>acc){
			//step++;
			dx=-B*grad;
			lambda=1.0;
			while(true){
				s=lambda*dx; 
				if(f(start+s)<f(start)){ // accept the step and update B
					start+=s;		
					gradOld=grad;
					grad=gradient(f,start);
					y=grad-gradOld;
					u=s-B*y;
					uy=u.dot(y);
					if(Abs(uy)>MEPS){
						dB=matrix.outer(u,u)/uy;
						B+=dB;
					}
					break;
				}
				lambda/=2;
				if(lambda < 1.0/1024){
					//WriteLine($"backtracking at {it}");
					start+=lambda*dx;
					grad=gradient(f,start);
					B.set_identity();
					break;
				}
			}
		}
		return start;
	}
	
	public static vector fit(Func<vector, double> f, vector guess, genlist<double> xs, genlist<double> ys, genlist<double> errs, double acc){ // vectros will not work, have to use List<double> for ex yx and errs
		Func<vector, double> dev = (x) => { // dev = deviation function
			double sum = 0; for(int i=0; i<xs.size; i++) sum += Pow((f(new vector(xs[i], x[0], x[1], x[2])) - ys[i])/errs[i], 2); return sum;
		};
		return qnewton(dev, guess, acc);
	}

	public static vector gradient(Func<vector, double> f, vector x){
		int n = x.size;
		double dx;
		vector xnew=x.copy();
		vector grad = new vector(n);
		for(int i=0;i<n;i++){
			dx=Abs(x[i])*Pow(2,-26);
			xnew[i]+=dx;
			grad[i]=(f(xnew)-f(x))/dx;
			xnew[i]-=dx;
		}
		return grad;
	}
	

	static public double ros(vector v){
		it++;
		double res=Pow(1-v[0],2)+100*Pow(v[1]-Pow(v[0],2),2);
		return res;
	}
	static public double him(vector v){
		it++;
		double x = v[0];
		double y = v[1];
		double res=Pow(Pow(x,2)+y-11,2)+Pow(x+Pow(y,2)-7,2);
		return res;
	}

	static public double test1(vector v){
		it++;
		double res=Pow(v[0],2)+Pow(v[1],2);
		return res;
	}

	static public double test2(vector v){
		it++;
		double res=Pow(v[0]+4,2)+Pow(v[1]+4,2);
		return res;
	}
		
	static public double test3(vector v){
		it++;
		double res=v[0]+Pow(v[1],2);
		return res;
	}

	static public double test4(vector v){ // f(x, y) = (x - 1)^4 + (x - 1)^2 + y^2
		it++;
		double x = v[0];
		double y = v[1];
		double res = Pow(x-1,4)+Pow(x-1,2)+Pow(y,2); // local minima at (1, 0) and (-1, 0)
		return res;
	}


	static public double test5(vector x){ // f(x, y) = x^4-7*x^3-32*x^2+y^2
		it++;
		return Pow(x[0],4)-7*Pow(x[0],3)-32*Pow(x[0],2)+Pow(x[1],2);
	}




	public static void Main(){
		// A
		WriteLine("==============================[A]==============================");
		vector start;
		WriteLine("test of x^2+y^2");
		start=new vector(10,-10);
		qnewton(test1,start).print("awnser: ");
		WriteLine($"should be approx (0,0) the function was called {it} times\n"); it=0;
		
		WriteLine("test of (x+4)^2+(y+4)^2");
		start=new vector(5,-7);
		qnewton(test2,start).print("awnser: ");
		WriteLine($"should be approx (-4,-4) the function was called {it} times\n"); it=0;

		WriteLine("test of f(x, y) = x^4-7*x^3-32*x^2+y^2");
		start=new vector(10,-2);
		qnewton(test5,start).print("awnser: ");
		WriteLine($"should be approx (7.41,0) the function was called {it} times\n"); it=0;

		start = new vector(3,3);
		WriteLine($"minimum of the Rosenbrock's valley function");
		qnewton(ros,start).print("awnser: ");
		WriteLine($"should be approx (1,1) the function was called {it} times\n"); it=0;
		WriteLine();
				
		WriteLine("minimum of the Himmelblau's function, has four local minima at (3,2), (-2.8,3.1) (-3.8,-3.3) and (3.6,-1.8):\n");
		qnewton(him,new vector( 2.5  ,  2.5)).print("awnser: ");
		WriteLine($"should be approx (3,2) the function was called {it} times\n"); it=0;
		qnewton(him,new vector(-2.5  ,  2.8)).print("awnser");
		WriteLine($"should be approx (-2.8,3.1) the function was called {it} times\n"); it=0;
		qnewton(him,new vector(-3.5  , -3.0)).print("awnser: ");
		WriteLine($"should be approx (-3.8,3.1) the function was called {it} times\n"); it=0;
		qnewton(him,new vector( 3.5  , -1.3)).print("awnser: ");
		WriteLine($"should be approx (3.6,-1.9) the function was called {it} times\n"); it=0;
		
		// B
		WriteLine("==============================[B]==============================");
		Func<vector, double> bw = (x) => x[1]/(Pow(x[0]-x[2], 2) + Pow(x[3],2)/4); // E = 1, m = 2, gamma = 3, A = 4
		var energy = new genlist<double>();
		var signal = new genlist<double>();
		var error  = new genlist<double>();
		
		string[] data = File.ReadAllLines("higgs.data");
		string[] words;
		var separators = new char[] {' ','\t'};
		var options = StringSplitOptions.RemoveEmptyEntries;
		for(int i = 0; i<data.Length; i++){
			words = data[i].Split(separators, options);
			energy.add(double.Parse(words[0]));
			signal.add(double.Parse(words[1]));
			error .add(double.Parse(words[2]));
		}
		
		WriteLine("energy\tsignal\terror");
		for(int i=0;i<energy.size;i++){
			WriteLine($"{energy[i]}\t{signal[i]}\t{error[i]}");
		}
	
		vector guess = new vector(2,120,2);
		vector bw_params = fit(bw, guess, energy, signal, error,1e-4);
		

		double A,m,gamma;A=bw_params[0];m=bw_params[1];gamma=bw_params[2];
		WriteLine($"parameters for the Breit-Wigner function: A={A}\tm={m}\tgamma={gamma}");
				
		double bwData;
		string toWrite="";
		for(double E=100;E<160;E+=0.1){
			bwData = bw(new vector(E,A,m,gamma));
			toWrite+=$"{E}\t{bwData}\n";
			
		}
		File.WriteAllText("bwfit.data",toWrite);
		
		/*
		qnewton(him,new vector(-2.8, 3.1)).print($"should be approx (-2.8,3.1):"); 
		WriteLine($"it took {step} steps"); step=0;
		qnewton(him,new vector(-3.8,-3.3)).print($"should be approx (-3.8,-3.3):"); 
 		WriteLine($"it took {step} steps"); step=0;
		qnewton(him,new vector( 3.6,-1,8)).print($"should be approx (3.6,-1.8):"); 
		WriteLine($"it took {step} steps"); step=0;
		*/
		/* // testus maximus
		it=0;
		qnewton(test5, new vector(10,-2),0.001).print("test5");
		WriteLine($"it took {it} steps"); it=0;
		*/
	}	
}// class


using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class main{
	
	static double f1(double t){
		return 1.0;
	}

	static double f2(double t){
		return t;
	}

	static (vector, matrix) lsfit(Func<double, double>[] fs, vector x, vector y, vector dy){
		int k = fs.Length;
		int i = x.size;
		vector b = new vector(i);
		matrix A = new matrix(i,k);
		for(int n=0; n<i; n++){
			b[n] = y[n]/dy[n];
			for(int m=0; m<k; m++){
				A[n,m]=fs[m](x[n])/dy[n];
			}
		}
		var R = new matrix(k,k);
		var Q = A.copy();
		QRGS.decomp(Q,R);
		vector c = QRGS.solve(Q,R,b);
		matrix invR = QRGS.inverseR(R);
		matrix S = invR*invR.transpose();
		return (c,S);


		
	}


	public static void Main(){
		double[] ts 		= {1  ,2  ,3 ,4 ,6 ,9  ,10   ,13  ,15};
		vector t = new vector(ts);
		
		double[] acts 		= {117,100,88,72,53,29.5,25.2,15.2,11.1};
		vector act = new vector(acts);
		vector logAct = act.copy();
		for(int i=0; i<act.size; i++){
			logAct[i]=Log(act[i]);
		}
		
		double[] errActs 	= {5  ,5  ,5 ,4 ,4 ,3   ,3   ,2   ,2};
		vector errAct = new vector(errActs);
		vector logErrAct = errAct.copy();
		for(int i=0; i<act.size; i++){
			logErrAct[i]=errAct[i]/act[i];
		}
		
		//make data
		string toWrite="";
		for(int i=0; i<act.size; i++){
			logErrAct[i]=errAct[i]/act[i];
			toWrite+=$"{t[i]}\t{act[i]}\t{errAct[i]}\n";
		}
		File.WriteAllText("radium.data",toWrite);	
	
		
		WriteLine("====================[A]====================");
		WriteLine("c1 and c2");
		Func<double, double>[] fs = {f1,f2};
		var solution = lsfit(fs,t,logAct,logErrAct).Item1;
		solution.print();
		double halflife = -Log(2)/solution[1];
		WriteLine($"halflife: {halflife} [days]");
		WriteLine("halflife of radium is 3.6 days");
		
		//make fit data
		toWrite="";
		double res = 100;
		for(int i=0;i<res;i++){
			double ti = i*1.0/res*t[t.size-1];
			double ft = Exp(solution[0])*Exp(solution[1]*ti);
			toWrite+=$"{ti}\t{ft}\n";
		}
		File.WriteAllText("radiumfit.data",toWrite);

		
		WriteLine("====================[B]====================");	
		var S = lsfit(fs,t,logAct,logErrAct).Item2;
		WriteLine("Covariance matrix");
		S.print();
		double uncer = Log(2)/Pow(solution[1],2)*Sqrt(S[1][1]); /* gange på hældningen */
		WriteLine($"halflife and its uncertainty in days");
		WriteLine($"{halflife} +-{uncer}");






		






	}	
}// class


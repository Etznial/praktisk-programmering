using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;
using System.Diagnostics;

class main{
	public static int it=0;
	public static int step=0;
	public static vector qnewton(
		Func<vector,double> f, 	/* objective function */
		vector start, 		/* starting point */
		double acc = 1e-3		/* accuracy goal, on exit |gradient| should be < acc */){
		
		int n = start.size;
		double lambda=1;
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
			if(Abs(start[i]) < Pow(2,-26)) start[i] = Pow(2,-26);
		}
		while(grad.norm() > acc){
			step++;
			dx=-B*grad;
			while(true){
				s=lambda*dx; 
				if(f(start+s)<f(start)){
					start+=s;		
					gradOld=grad;
					grad=gradient(f,start);
					y=grad-gradOld;
					u=s-B*y;
					uy=u.dot(y);
					//check if ny is zero
					if(Abs(uy) > 1e-6){
						dB=matrix.outer(u,u)/uy;
						B+=dB;
					}
					break;
				}
				lambda/=2;
				if(lambda < 1.0/1024){
					start+=s;
					grad=gradient(f,start);
					B.set_identity();
					break;
				}
			}
		}
		return start;
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

	public static vector qnewtonOld(
		Func<vector,double> f, 	/* objective function */
		vector start, 		/* starting point */
		double acc = 2e-23		/* accuracy goal, on exit |gradient| should be < acc */){
		int n = start.size;
		/* Hessian matrix */
		matrix H = new matrix(n,n);
		vector x = start.copy();
		vector x1 = start.copy();
		vector x12 = start.copy(); /* x1 = xdx1, x12 = xdx1dx2 */
		double dx1, dx2;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				dx1 = Abs(x1[i])*Pow(2,-26);
				dx2 = Abs(x12[j])*Pow(2,-26);
				x1[i]+=dx1; x12[i]+=dx1; x12[j]+=dx2;/* fix x'erne */
				WriteLine($"f(x12)={f(x12)}\tf(x1)={f(x1)}\tf(x)={f(x)}");
				x.print("x  =");
				x1.print("x1 =");
				x12.print("x12=");
				H[i,j]=((f(x12)-f(x))/dx1 - (f(x1)-f(x))/dx1)/dx2;
				x1[i]-=dx1; x12[i]-=dx1; x12[j]-=dx2;/* fix x'erne */
			}
		}
		H.print("Hessian matrix:");
		


		var temp = new vector(n);
		return temp;

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


	public static void Main(){
		var start = new vector(2,2);
		WriteLine($"minimum of the Rosenbrock's valley function, it took {step} steps, has one global minima at (1,1)");
		qnewton(ros,start).print("awnser: "); step=0;
		WriteLine();
				
		WriteLine("minimum of the Himmelblau's function, has four local minima at (3,2), (-2.8,3.1) (-3.8,-3.3) and (3.6,-1.8):\n");
		qnewton(him,new vector( 2.5  ,  2.5)).print("awnser: ");
		WriteLine($"should be approx (3,2) it took {step} steps\n"); step=0;
		qnewton(him,new vector(-2.5  ,  2.8)).print("awnser");
		WriteLine($"should be approx (-2.8,3.1) it took {step} steps\n"); step=0;
		qnewton(him,new vector(-3.5  , -3.0)).print("awnser: ");
		WriteLine($"should be approx (-3.8,3.1) it took {step} steps\n"); step=0;
		qnewton(him,new vector( 3.5  , -1.3)).print("awnser: ");
		WriteLine($"should be approx (3.6,-1.9) it took {step} steps\n"); step=0;




		/*
		qnewton(him,new vector(-2.8, 3.1)).print($"should be approx (-2.8,3.1):"); 
		WriteLine($"it took {step} steps"); step=0;
		qnewton(him,new vector(-3.8,-3.3)).print($"should be approx (-3.8,-3.3):"); 
 		WriteLine($"it took {step} steps"); step=0;
		qnewton(him,new vector( 3.6,-1,8)).print($"should be approx (3.6,-1.8):"); 
		WriteLine($"it took {step} steps"); step=0;
		*/
	}	
}// class


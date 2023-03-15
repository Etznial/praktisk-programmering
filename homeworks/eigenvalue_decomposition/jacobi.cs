using System;
using static System.Console;
using static System.Math;
public static class jacobi {
	public static void timesJ(matrix A, int p, int q, double theta){ /* A -> A*J */
		double c = Cos(theta), s=Sin(theta);
		for(int i=0; i<A.size1; i++){
			double aip=A[i,p],aiq=A[i,q];
			A[i,p]=c*aip-s*aiq;
			A[i,q]=s*aip+c*aiq;
		}
	}
	public static void Jtimes(matrix A, int p, int q, double theta){
		double c = Cos(theta), s=Sin(theta);
		for(int j=0; j<A.size1; j++){
			double apj=A[p,j], aqj=A[q,j];
			A[p,j]= c*apj+s*aqj;
			A[q,j]=-s*apj+c*aqj;
		}
	}
	public static void decomp(matrix D, matrix V){
		bool changed;
		int n = D.size1;
		V.set_identity();
		//int it=0;
		//int itch=0;
		do{
			changed = false;
			for(int p=0;p<n-1;p++){
				for(int q=p+1;q<n;q++){
					double apq=D[q,p], app=D[p,p], aqq=D[q,q];
					double theta = 0.5*Atan2(2*apq,aqq-app);
					if(apq > 1e-15 | -apq >1e-15){ // do rotation
						changed = true;
						jacobi.Jtimes(D,p,q,-theta); // D←JT*D
						jacobi.timesJ(D,p,q, theta); // D←D*J 
						jacobi.timesJ(V,p,q, theta); // V←V*J
						
						
						//D.print();
						//itch++;
					}
				}
			}
			//it++;
		}while(changed);
	}

} // jacobi



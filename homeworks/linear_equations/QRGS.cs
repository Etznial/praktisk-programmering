using System;
using static System.Console;
using static System.Math;

public class QRGS{
	public static void decomp(matrix Q, matrix R){ 
/* A must be square and invertable, A=QR where R should be an upper triangular matrix, Q should be and orthoganol or unitary matrix -> Q^T=Q^(-1). The columbs of Q form an orthonalmal basis*/
		int m = Q.size2; /* width of matrix */
		for(int i=0; i<m; i++){
			R[i][i]=Q[i].norm();
			Q[i]/=R[i,i];
			for(int j=i+1;j<m;j++){
				R[j][i]=Q[i].dot(Q[j]);
				Q[j]-=Q[i]*R[j][i];
			}
		}
	}
	public static vector solve(matrix Q, matrix R, vector b){ /* with back substitution, solver ting på formen Rx=c, c=QTb */
		int n = b.size;
		vector solutions = new vector(n);
		for(int i=n-1; i>=0;i--){
			double tempSum = 0;
			for(int k=i+1; k<n; k++){
				tempSum += R[k][i]*solutions[k];
			}
		solutions[i]=1.0/R[i][i]*(b[i]-tempSum);
		}
		return solutions; 
	}
		
	public static double det(matrix R){ /* R skal være decomposed inden man indsætter den */
		int m = R.size2;
		double tempSum = 1;
		for(int i=0; i<m; i++){
			tempSum = tempSum*R[i][i];
		}
		return tempSum;
	}

	public static matrix inverse(matrix Q, matrix R){
		int n = R.size1;
		matrix AI = new matrix(n,n); /* Detter er A inverse */
		for(int i=0; i<n; i++){
			vector e = new vector(n);
			e[i]=1;
			vector QTe = Q.transpose()*e;
			AI[i] = solve(Q,R,QTe);		
		}
		return AI;
	}

}


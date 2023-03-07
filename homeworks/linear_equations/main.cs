using System;
using static System.Console;
using static System.Math;

public static class QRGS{
	public static void decomp(matrix A, matrix R){ 
/* A must be square and invertable, A=QR where R should be an upper triangular matrix, Q should be and orthoganol or unitary matrix -> Q^T=Q^(-1). The columbs of Q form an orthonalmal basis*/
		int m = A.size1;
		matrix Q = A.copy();
		R = new matrix(m,m);
		for(int i=0; i<m; i++){
			R[i,i]=Q[i].norm();
			Q[i]/=R[i,i];
			for(int j=i+1;j<m;j++){
				R[i,j]=Q[i].dot(Q[j]);
				Q[j]-=Q[i]*R[i,j];
			}
		}
	}
}


class main{
	public static void Main(){
		WriteLine("====================[start]====================");
		WriteLine("create vector");
		var a = new vector(3);
		a.print();
		WriteLine("set values of vector");
		for(int i=0;i<a.size;i++){
			a[i]=i+1;	
		}
		a.print();
		WriteLine("create matrix");
		var A = new matrix(3,3);
		A.print();
		WriteLine("set values of matrix");
		int k = 1;
		for(int c=0; c<A.size1;c++){
			for(int r=0; r<A.size2;r++){
				A[r][c] = k;
				k++;					
			}
		}
		A.print();
		WriteLine("====================[end]====================");

	}
}// class


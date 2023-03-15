using System;
using static System.Console;
using static System.Math;
using static System.Diagnostics;

public static class QRGS{
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


class main{
	public static void Main(){
		var random = new System.Random(1);
		{// A
		WriteLine("====================[A]====================");
	
	
	
		
		
		
		
		{//decomp
		WriteLine("====================[decomp]====================");
		WriteLine("=====[generate a random number (n>m) matrix A (of a modest size)]=====");
		int n = random.Next(7,9);
		int m = random.Next(4,6);
		WriteLine($"n = {n}, m = {m}");
		//WriteLine("generate random matrix");
		matrix A = new matrix(n,m);
		
                //WriteLine("set matrix values with random doubles");
	       	for(int c=0; c<A.size2;c++){
			for(int r=0; r<A.size1;r++){
				A[c][r] = random.NextDouble();
			}
		}	
		WriteLine("Print matrix A");
		A.print();
		WriteLine("=====[factorize matrix A it into QR]=====");
		var R = new matrix(m,m); /* making R with dimensions [m,m] */
		matrix Q = A.copy(); /* making Q as a copy of A*/
		QRGS.decomp(Q,R); 
		WriteLine("Print matrix Q");
		Q.print(); 
		WriteLine("Print matrix R");
		R.print();
		WriteLine("=====[check that R is upper triangular]=====");
		WriteLine("See the print of R just above");
		WriteLine("=====[check that QTQ=1]=====");
		WriteLine();
		matrix QTQ = Q.transpose()*Q;
		WriteLine("Print  QTQ = Q.transpose()*Q");
		QTQ.print();
		WriteLine("=====[check that QR=A]=====");
		matrix QR = Q*R;
		WriteLine("Print of matrix Q*R");
	       	QR.print();	
		WriteLine("Print of matrix A");
		A.print();
		}//decomp
		/*
		WriteLine("====================[start of test of det]====================");
		WriteLine("detA = detQ*detR = +-detR = produktet af diagonalen af R");
		WriteLine($"detR = {QRGS.det(R)}");
		WriteLine("====================[end of test of det]====================");
		*/
		{ // solve
		WriteLine("====================[solve]====================");
		WriteLine("=====[generate a random square matrix A (of a modest size)]=====");
		int n = random.Next(7,9);
		int m = n;
		matrix A = new matrix(m,n);
		for(int c=0; c<n;c++){
			for(int r=0; r<m;r++){
				A[c][r] = random.NextDouble();
			}
		}
		WriteLine("Print matrix A");
		A.print();
		
		WriteLine("=====[generate a random vector b (of the same size)]=====");
		vector b = new vector(m);
		for(int i=0;i<b.size;i++){
			b[i]=random.NextDouble();	
		}
		WriteLine("Print vector b");
		b.print();
		WriteLine("=====[factorize matrix A it into QR]=====");
		var R = new matrix(n,n); /* making R with dimensions [m,m] */
		matrix Q = A.copy(); /* making Q as a copy of A*/
		QRGS.decomp(Q,R); 
		WriteLine("Print matrix Q");
		Q.print(); 
		WriteLine("Print matrix R");
		R.print();

		WriteLine("=====[solve QRx=b]=====");
		vector QTb = Q.transpose()*b;
		WriteLine("Print vector QTb");
		QTb.print();
		WriteLine("Solutions for x:");
		vector xs = QRGS.solve(Q,R,QTb);
		xs.print();
		


		WriteLine("=====[check that Ax=b]=====");
		WriteLine("Print vector b");
		b.print();
		WriteLine("Print vector Ax");
		vector Ax = A*xs;
		Ax.print();
		}// solve
	        }// A
		{// B
		WriteLine("====================[B]====================");
		WriteLine("=====[generate a random square matrix A (of a modest size)]=====");
		int n = random.Next(7,9);
		int m = n;
		matrix A = new matrix(m,n);
		for(int c=0; c<n;c++){
			for(int r=0; r<m;r++){
				A[c][r] = random.NextDouble();
			}
		}
		WriteLine("Print matrix A");
		A.print();
		WriteLine("=====[factorize matrix A it into QR]=====");
		var R = new matrix(n,n); /* making R with dimensions [m,m] */
		matrix Q = A.copy(); /* making Q as a copy of A*/
		QRGS.decomp(Q,R); 
		WriteLine("Print matrix Q");
		Q.print(); 
		WriteLine("Print matrix R");
		R.print();
		WriteLine("=====[calculate the inverse B]=====");
		matrix B = QRGS.inverse(Q,R);
		WriteLine("Print matrix B");
		B.print();
		WriteLine("=====[check that AB=I, where I is the identity matrix]=====");
		WriteLine("Print matrix A*B");
		matrix I = A*B;
		I.print();
		}// B
		{// C
		WriteLine("====================[C]====================");
		
		Stopwatch stopwatch = new Stopwatch();

		for(int i=0; i<N; i++){
			stopwatch.Start();

			make NxN matrix
			
			give random numbers
					
		    	fac NxN matrix 
			stopwatch.Stop();
			TimeSpan elapsed = stopwatch.Elapsed;
		}
		

		
		WriteLine($"{N} {time.elapsed}");
		

		}// C
	}
}// class


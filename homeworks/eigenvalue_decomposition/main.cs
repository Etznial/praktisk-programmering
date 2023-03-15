using System;
using static System.Console;
using static System.Math;
class main{
	public static void Main(string[] args){
		double rmax = 0;
		int N = 0;

		foreach(var arg in args){
			var words = arg.Split(':');
			if(words[0] == "-rmax"){rmax = double.Parse(words[1]);}
			if(words[0] == "-N"){N = int.Parse(words[1]);}
		}
		double dr = rmax/N;
		WriteLine($"rmax:{rmax} N:{N} dr:{dr}");
		/* create hamiltonian */
		int npoints = (int)(rmax/dr)-1;
		vector r = new vector(npoints);
		for(int i=0; i<npoints; i++)r[i]=dr*(i+1);
		matrix H = new matrix(npoints,npoints);
		for(int i=0; i<npoints-1; i++){
			H[i,i]	=-2;
			H[i,i+1]= 1;
			H[i+1,i]= 1;
		}
		H[npoints-1,npoints-1]=-2;
		matrix.scale(H,-0.5/(dr*dr));
		for(int i=0; i<npoints; i++)H[i,i]+=-1/r[i];
		
		/* create symetric matrix */
		var random = new System.Random(1);
		int n = random.Next(5,5);
		matrix A = new matrix(n,n);
		for(int i=0; i<A.size1; i++){
			for(int j=A.size1-1; j>=i; j--){
				A[i,i]=random.NextDouble();
				A[i,j]=random.NextDouble();
				A[j,i]=A[i,j];
			}
		}
		matrix D = A.copy();
		matrix V = new matrix(n,n);
		matrix I = V.copy();
		I.set_identity();
		jacobi.decomp(D,V);

		//A.print();
	
		WriteLine("==============================[A]==============================");
		WriteLine("print A");
		A.print();
		WriteLine("print D");
		D.print();
		WriteLine("print V");
		V.print();

		
		matrix VTAV = V.transpose()*A*V;
		WriteLine($"print VT*A*V==D is {D.approx(VTAV)}");
		WriteLine("print VTAV");
		VTAV.print();
		
		matrix VDVT = V*D*V.transpose();
		WriteLine($"V*D*VT==A is {A.approx(VDVT)}");
		WriteLine("print VDVT");
		VDVT.print();
		
		matrix VTV = V.transpose()*V;
		WriteLine($"VT*V==1 is {VTV.approx(I)}");
		WriteLine("print VTV");
		VTV.print();
		
		matrix VVT = V*V.transpose();
		WriteLine($"V*VT==1 is {VVT.approx(I)}");
		WriteLine("print VVT");
		VVT.print();
	
		WriteLine("==============================[B]==============================");
		WriteLine("print H");
		H.print();	
		D = H.copy();
		n = H.size1;
		V = new matrix(n,n);
		jacobi.decomp(D,V);
		WriteLine("print D");
		D.print();
		
		double lowest = D[0,0];
		for(int i=0; i>n-1; i++){
			if(lowest>D[i,i]){
				lowest=D[i,i];
			}
		}
		WriteLine($"lowest eigen-function:{lowest}");
		WriteLine("print V");
		V.print();
		
	}
}// class


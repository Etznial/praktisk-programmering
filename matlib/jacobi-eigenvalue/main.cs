using System;

public class test{

public static void Main(string[] args){
int n=5, max_print=8;
if (args.Length > 0) n = int.Parse(args[0]);
Console.WriteLine("n={0}",n);
int nvals=n;
if(args.Length>1) nvals = int.Parse(args[1]);

var rnd = new Random(1);
var a = new matrix(n,n);
vector e = new vector(n);
matrix v = new matrix(n,n);
if(n>max_print)
	{
	for(int i=0;i<n;i++)for(int j=i;j<n;j++)
		a[i,j]=a[j,i]=2*(rnd.NextDouble()-0.59);
	int r=jacobi.cyclic(a,e,nvals,V:v);
	Console.WriteLine("sweeps={0} e[0]={1}",r,e[0]);
	return;
	}
else
	{
	for(int i=0;i<n;i++)for(int j=i;j<n;j++)
		{
		a[i,j]=a[j,i]=rnd.NextDouble();
		}
	Console.WriteLine("Eigenvalue Decomposition A=V*D*V^T");
	a.print("Random matrix A:");
	matrix b = a.copy();
	int sweeps=jacobi.cyclic(a,e,nvals,v);
	System.Console.WriteLine("Number of sweeps={0}",sweeps);
	var tmp=(v.T*b*v);
	tmp.print("V^T*A*V (should be diagonal):");
	e.print("Eigenvalues (should equal the diagonal elements above):\n");
	}
}
}

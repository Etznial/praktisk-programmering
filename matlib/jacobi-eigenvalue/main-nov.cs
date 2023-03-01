using System;

public class test{

public static void Main(string[] args){
int n=5;
if (args.Length > 0) n = int.Parse(args[0]);
Console.WriteLine("n={0}",n);

var rnd = new Random(1);
var a = new matrix(n,n);
vector e = new vector(n);
	for(int i=0;i<n;i++)for(int j=i;j<n;j++)
		a[i,j]=a[j,i]=2*(rnd.NextDouble()-0.59);
	int r=jacobi.cyclic(a,e,V:null);
	Console.WriteLine("sweeps={0} e[n-1]={1}",r,e[n-1]);
	return;
}
}

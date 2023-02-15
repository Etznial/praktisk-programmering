using System;
using static System.Console;
using static System.Math;
public static class main{
	public static void Main(){
		int n=9;
		double[] a = new double[n];
		for(int i=0;i<n;i++) Write($"a[{i}]={a[i]} ");
		WriteLine();
		for(int i=0;i<n;i++) a[i]=i; 
		for(int i=0;i<n;i++) Write($"a[{i}]={a[i]} ");
		WriteLine();
		WriteLine($"array length = {a.Length}");
		foreach(double ai in a) Write($"ai = {ai} ");
		WriteLine();
	
	}
}

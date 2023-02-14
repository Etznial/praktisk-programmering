using System;
using static System.Console;
using static System.Math;

public static class main{
	public static void Main(){
		
		/*1*/
		
		int i=1; while(i+1>i) {i++;}
		WriteLine($"My max int = {i}");
		WriteLine($"MaxValue = {int.MaxValue}");

		int j=-1; while(j-1<j) {j--;}
		WriteLine($"My min int = {j}");
		WriteLine($"MinValue = {int.MinValue}");
		
		/*2*/

		double x=1; while(1+x!=1){x/=2;} x*=2;
		WriteLine($"Mindste double computeren kan adderer 1 med uden at få 1 tilbage {x}");
		WriteLine($"Machine epsilon double = {Pow(2,-52)}");
		

		float y=1F; while (1+y!=1){y/=2;} y*=2;
		WriteLine($"Mindste float computeren kan adderer 1 med uden at få 1 tilbage {y}");
		WriteLine($"Single persision float = {Pow(2,-23)}");

		/*3*/ 
		
		int n = (int)1e6;
		double epsilon = Pow(2,-52);
		double tiny = epsilon/2;
		double sumA = 0, sumB = 0;

		sumA+=1; for(int k=0;i<n;k++){sumA+=tiny;}
		for(int k=0;k<n;k++){sumB+=tiny;} sumB+=1;

		WriteLine($"sumA-1 = {sumA-1:e} should be {n*tiny:e}");
		WriteLine($"sumB-1 = {sumB-1:e} should be {n*tiny:e}");
		WriteLine("Addere man noget lille til 1 for man 1 tilbage, men adderer man først en masse tiny sammen så er tiny stor nok til at kunne adderes med 1 uden kun at få 1 tilbage");
	
		/*4*/

		double d1 = 0.1+0.1+0.1+0.1+0.1+0.1+0.1+0.1;
		double d2 = 8*0.1;

		WriteLine($"d1={d1:e15}");
		WriteLine($"d2={d2:e15}");
		WriteLine($"d1==d2 ? => {d1==d2}");

		WriteLine($"sfuns.approx(d1,d2) ? => {sfuns.approx(d1,d2)}");
		
	
	
	}
}

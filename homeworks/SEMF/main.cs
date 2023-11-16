using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class main{
	// s. 261
	static double a_v=15.8;
	static double a_s=17.8;
	static double a_c=0.71;
	static double a_a=94.8;
	static double a_p=11.2;
	static double m_p=938; // MeV
	static double m_n=940; // MeV
	static double m_e=0.5; // MeV
	public static double f1(double A){
		return a_v*A;
	}
	public static double f2(double A){
		return -a_s*Pow(A,2.0/3.0);
	}
	public static double f3(double A, double Z){
		return -a_c*Z*(Z-1)*Pow(A,-1.0/3.0);
	}
	public static double f4(double A, double Z){
		return -a_a*Pow(Z-A/2.0,2)/A;
	}
	public static double f5(double A, double Z){
		if(A%2==0 && Z%2==0){
			return a_p*Pow(A,-1.0/2.0);
		}
		if(A%2!=0 && Z%2!=0){
			return -a_p*Pow(A,-1.0/2.0);
		}
		else{
			return 0;
		}
	}
	public static double B(double Z,double A){
		return f1(A)+f2(A)+f3(A,Z)+f4(A,Z)+f5(A,Z);
	}
	public static double m_A(double Z,double A){
		return Z*(m_p+m_e)+(A-Z)*m_n-B(Z,A);
	}
	
	public static void Main(){
		WriteLine("\tBinding Energy [MeV]\t\tBinding Energy [keV] IAEA");
		// 37Cl A = 37-17 = 20, Z = 17
		WriteLine($"37Cl:\t{B(37,20)}\t\t{37*8570.2816}");
		// 56Fe A = 56-26 = 30, Z = 26
		WriteLine($"56Fe:\t{B(56,30)}\t\t{56*8790.3563}");

		//make data for 37Li
		string toWrite="";
		int A=37;
		for(int i=0;i<A;i++){
			toWrite+=$"{i}\t{B(i,A)}\n";
		}
		File.WriteAllText("BE37Cl.data",toWrite);
		
		//make data SEMF
		toWrite="";
		for(int i=0;i<A;i++){
			toWrite+=$"{i}\t{m_A(i,A)}\n";
		}
		File.WriteAllText("ME37Cl.data",toWrite);
			
	}
// class
}

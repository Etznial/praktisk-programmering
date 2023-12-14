using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class main{
	// s. 261
	static double a_v = 15.56;
	static double a_s = 17.23;
	static double a_c = 0.697;
	static double a_a = 93.14;
	static double a_p = 12.0;
	static double m_n = 939.56542052;
	static double m_p = 938.27208816;
	static double m_e = 0.51099895000;
	public static double f1(double A){
		return a_v*A;
	}
	public static double f2(double A){
		return -a_s*Pow(A,2.0/3.0);
	}
	public static double f3(double A, double Z){
		return -a_c*Z*(Z-1.0)*Pow(A,-1.0/3.0);
	}
	public static double f4(double A, double Z){
		return -a_a*Pow(Z-A/2.0,2.0)/A;
	}
	public static double f5(double A, double Z){
		if(A%2==0 && Z%2==0){
			return a_p*Pow(A,-1.0/2.0);
		}
		if(A%2!=0 && Z%2!=0){
			return -a_p*Pow(A,-1.0/2.0);
		}
		else{
			return 0.0;
		}
	}
	public static double B(double Z,double A){
		return f1(A)+f2(A)+f3(A,Z)+f4(A,Z)+f5(A,Z);
	}
	public static double SEMF(double Z,double A){
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
			toWrite+=$"{i}\t{SEMF(i,A)}\n";
		}
		File.WriteAllText("ME37Cl.data",toWrite);
			
	}
// class
}

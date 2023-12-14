using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class semf{
	static double mn = 939.56542052;
	static double mp = 938.27208816;
	static double me = 0.51099895000;

	static double a_v = 15.56;
	static double a_s = 17.23;
	static double a_c = 0.697;
	static double a_a = 93.14;
	static double a_p = 12.0;
	

	public static double alpha(double A){
		return mn - a_v + Pow(a_s/A,(1.0/3.0)) + a_a/4.0;
	}
	public static double beta(){
		return a_a + mn - mp - me;
	}
	public static double gamma(double A){
		return a_a/A + a_c/Pow(A,(1.0/3.0));
	}
	public static double delta(double A, double Z){
		if(A%2==0 && Z%2==0){return -a_p;}
		if(A%2!=0 && Z%2!=0){return  a_p;}
		else{return 0.0;}
	}
	public static double SEMF(double A, double Z){ 
		return alpha(A)*A-beta()*Z+gamma(A)*Pow(Z,2)+delta(A,Z)/Pow(A,(1/2));
	}
	//print("235U -> 87Br + 145La + 3n. Q-value:", SEMF(235, 92) - (SEMF(87, 35) + SEMF(145, 57) + 3*mn), "MeV")
	public static void Main(){
		//make data SEMF
		int A=37;
		string toWrite="";
		for(int i=0;i<A;i++){
			toWrite+=$"{i}\t{SEMF(A,i)}\n";
		}
		File.WriteAllText("ME37Clpy.data",toWrite);
		
	}
}

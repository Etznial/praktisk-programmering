using System;
using static System.Console;
using static System.Math;
using static matrix;
using System.IO;

class main{
	public static void latticeData(double a, double b, double theta, string dataName, int itt, double offsetx=0, double offsety=0){
		string toWrite="";
		for(int i=0;i<itt;i++){
			for(int j=0;j<itt;j++){
				if(j%2==0){
				 	toWrite+=$"{i*a+offsetx}\t{j*b*Math.Sin(theta)+offsety}\n";
				}else{
				 	toWrite+=$"{i*a+b*Math.Cos(theta)+offsetx}\t{j*b*Math.Sin(theta)+offsety}\n";
				}
			}
		}
		File.WriteAllText(dataName,toWrite);
	}
	
	public static void Main(){

		//Make data
		double Au=0.288; // [Å]
		double NbS2=0.33; // [Å]
		latticeData(Au,Au ,PI/3,"testLattice1.data",50		,1		,1);
		latticeData(NbS2,NbS2,PI/3,"testLattice2.data",43	,3		,3);
		latticeData(NbS2,NbS2,PI/3,"testLattice3.data",43	,3+19*NbS2*0.5	,3+19*Tan(PI/6.0)*NbS2*0.5);
	}
// class
}

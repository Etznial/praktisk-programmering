using System;
using static System.Console;
using static System.Math;
class main{
	public static void Main(string[] args){
		complex minusEn = new complex(-1,0);
		WriteLine($"cmath.sqrt(-1) = {cmath.sqrt(minusEn)}");
		WriteLine($"Awnser from Wolframalpha: i");
		WriteLine($"cmath.sqrt(-1).approx(new complex(0,1)) = {cmath.sqrt(minusEn).approx(new complex(0,1))}\n");
		
		WriteLine($"cmath.sqrt(i) = {cmath.sqrt(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: 1/sqrt(2) + i/sqrt(2)");
		WriteLine($"cmath.sqrt(i).approx(new complex(1/cmath.sqrt(2), i/cmath.sqrt(2))) = {cmath.sqrt(cmath.I).approx(new complex(1/cmath.sqrt(2), 1/cmath.sqrt(2)))}\n");
		
		WriteLine($"cmath.exp(cmath.I) {cmath.exp(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: cos(1) + sin(1)i");
		WriteLine($"cmath.exp(cmath.I).approx(new complex(Cos(1), Sin(1))) = {cmath.exp(cmath.I).approx(new complex(Cos(1), Sin(1)))}\n");
		
		WriteLine($"cmath.exp(cmath.I*PI) {cmath.exp(cmath.I*PI)}");
		WriteLine($"Awnser from Wolframalpha: -1");
		WriteLine($"cmath.exp(cmath.I*PI).approx(-1) = {cmath.exp(cmath.I*PI).approx(-1)}\n");

		WriteLine($"cmath.pow(cmath.I,cmath.I) = {cmath.pow(cmath.I,cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: Exp(-PI/2)");
		WriteLine($"cmath.pow(cmath.I,cmath.I).approx(new complex(Exp(-PI/2),0)) = {cmath.pow(cmath.I,cmath.I).approx(new complex(Exp(-PI/2),0))}\n");
		
		WriteLine($"cmath.log(cmath.I) = {cmath.log(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: I*PI/2");
		WriteLine($"cmath.log(cmath.I).approx(new complex(0,PI/2)) = {cmath.log(cmath.I).approx(new complex(0,PI/2))}\n");
		
		WriteLine($"cmath.sin(cmath.I*PI) = {cmath.sin(cmath.I*PI)}");
		WriteLine($"Awnser calculated: i/2*(Exp(PI)-Exp(-PI))");
		WriteLine($"cmath.sin(cmath.I*PI).approx(new complex(0,1.0/2*(Exp(PI)-Exp(-PI)))) = {cmath.sin(cmath.I*PI).approx(new complex(0,1.0/2*(Exp(PI)-Exp(-PI))))}");
	}
}

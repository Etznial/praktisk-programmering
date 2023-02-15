using System;
using static System.Console;
using static System.Math;
class main{
	public static void Main(string[] args){
		complex minusEn = new complex(-1,0);
		WriteLine($"cmath.sqrt(-1) = {cmath.sqrt(minusEn)}");
		WriteLine($"Awnser from Wolframalpha: i");
		WriteLine($"cmath.sqrt(i) = {cmath.sqrt(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: 0.707 + 0.707i");
		WriteLine($"cmath.exp(cmath.I) {cmath.exp(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: 0.540 + 0.841i");
		WriteLine($"cmath.exp(cmath.I*PI) {cmath.exp(cmath.I*PI)}");
		WriteLine($"Awnser from Wolframalpha: -1");
		WriteLine($"{cmath.log(cmath.I)}");
		WriteLine($"Awnser from Wolframalpha: 1.57i");
		WriteLine($"{cmath.sin(cmath.I*PI)}");
		WriteLine($"Awnser from Wolframalpha: 11.5i");
	}
}

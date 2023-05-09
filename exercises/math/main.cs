using System;
using static System.Console;
using static System.Math;

public static class main{
	public static void Main(){
		double sqrt2 = Sqrt(2);
		WriteLine($"sqrt(2)^2 = {sqrt2*sqrt2} should equal 2");
		
		double pow15 = Pow(2,1.0/5);
		WriteLine($"2^1.5^5 = {Pow(pow15,5)} should equal 2");		

		double ePI = Pow(Math.E, Math.PI);
		WriteLine($"Log(e^PI) = {Log(ePI)} should equal PI");

                double PIe = Pow(Math.PI, Math.E);
                WriteLine($"Log(PI^e) = {Log(PIe)} should equal e");

		double gamma1 = sfuns.gamma(1);
		WriteLine($"gamma1 = {gamma1} should equal 1");

		double gamma2 = sfuns.gamma(2);
		WriteLine($"gamma2 = {gamma2} should equal 1");

		double gamma3 = sfuns.gamma(3);
		WriteLine($"gamma3 = {gamma3} should equal 2");

		double gamma10 = sfuns.gamma(10);
		WriteLine($"gamma10 = {gamma10} should equal 362880");

                double lngamma1 = sfuns.lngamma(1);
                WriteLine($"lngamma1 = {lngamma1} should equal 0");

                double lngamma2 = sfuns.lngamma(2);
                WriteLine($"lngamma2 = {lngamma2} should equal 0");

                double lngamma3 = sfuns.lngamma(3);
                WriteLine($"lngamma3 = {lngamma3} should equal 0.693147");

                double lngamma10 = sfuns.lngamma(10);
                WriteLine($"lngamma10 = {lngamma10} should equal 12.8018");
	}
}

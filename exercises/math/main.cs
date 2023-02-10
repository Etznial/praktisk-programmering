using System;
using static System.Console;
using static System.Math;

public static class main{
	public static void Main(){
		double sqrtTest = Sqrt(2);
		WriteLine($"sqrtTest = {sqrtTest}");
		
		double powTest = Pow(2,1.0/5);
		WriteLine($"powTest = {powTest}");		

		double eTest = Pow(Math.E, Math.PI);
		WriteLine($"eTest = {eTest}");

                double piTest = Pow(Math.PI, Math.E);
                WriteLine($"piTest = {piTest}");

		double gamma1 = sfuns.gamma(1);
		WriteLine($"gamma1 = {gamma1}");

		double gamma2 = sfuns.gamma(2);
		WriteLine($"gamma2 = {gamma2}");

		double gamma3 = sfuns.gamma(2);
		WriteLine($"gamma3 = {gamma3}");

		double gamma10 = sfuns.gamma(10);
		WriteLine($"gamma10 = {gamma10}");

                double lngamma1 = sfuns.lngamma(1);
                WriteLine($"lngamma1 = {lngamma1}");

                double lngamma2 = sfuns.lngamma(2);
                WriteLine($"lngamma2 = {lngamma2}");

                double lngamma3 = sfuns.lngamma(2);
                WriteLine($"lngamma3 = {lngamma3}");

                double lngamma10 = sfuns.lngamma(10);
                WriteLine($"lngamma10 = {lngamma10}");
	}
}

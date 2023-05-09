using System;
using static System.Console;
using static System.Math;
public static class main{
	public static void print(this double x, string s){ /* x.print("x="); */
		Write(s);WriteLine(x);
		}	
	public static void print(this double x){ /* x.print() */
		x.print("");
	}

	public static void Main(){
		vec u = new vec(1,2,3);
		vec u_clone = new vec(1,2,3);
		vec v = new vec(2,3,4);
		u.print("u = ");
		v.print("v = ");
		(u+v).print("u+v = ");
		(2*u).print("2*u = ");
		WriteLine($"u%v = {u%v}");
		WriteLine($"u.dot(v) = {u.dot(v)}");
		
		double x = 1.23;
		x.print("x = ");
		(x+5).print("x+5 = ");
		(5*x).print("5*x = ");
		WriteLine($"u.cross(v) = {u.cross(v)}");
		WriteLine($"u.norm() = {u.norm()}");
		/* testing the approx method */
		WriteLine($"u.approx(u_clone) = {u.approx(u_clone)}");
		WriteLine($"u.approx(v) = {u.approx(v)}");

	}
}

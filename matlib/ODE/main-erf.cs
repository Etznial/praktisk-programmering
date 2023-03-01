using System;
using static System.Console;
using static System.Math;
using System.Collections;
using System.Collections.Generic;

class main{

static bool approx(double a,double b, double acc=1e-6, double eps=1e-6){
	if( Abs(a-b)<acc )return true;
	if( Abs(a-b)<eps*Max(Abs(a),Abs(b)) ) return true;
	return false;
}

static void Main(){

	double a=-4,b=4;
	for(double z=a;z<=b;z+=1.0/8)
		WriteLine($"{z} {matlib.erf(z)}");

	WriteLine("\n");
	WriteLine("-2   -0.995322265");
	WriteLine("-1   -0.842700793");
	WriteLine("-0.5	-0.520499878");
	WriteLine("-0.2	-0.222702589");
	WriteLine(" 0.2	 0.222702589");
	WriteLine(" 0.5	 0.520499878");
	WriteLine(" 1    0.842700793");
	WriteLine(" 2    0.995322265");
}
}

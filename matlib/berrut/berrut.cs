using System;
using System.Diagnostics;
public static partial class matlib{

public static double berrut(double[] x, double[] y, double z){
	int sign=1;
	double p=0,q=0;
	for(int i=0;i<x.Length;i++){
		if(z==x[i])return y[i];
		p+=sign*y[i]/(z-x[i]);
		q+=sign/(z-x[i]);
		sign*=-1;
		}
	return p/q;
}


}//matlib

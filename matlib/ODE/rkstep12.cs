using System;
public static partial class ODE{

public static (vector,vector) rkstep12
(Func<double,vector,vector> F, double x, vector y0, double h)
{// Runge-Kutta mid-point step
	vector k0  = F(x,y0);
	vector k12 = F(x+h/2,y0+k0*(h/2));
	vector yh=y0+k12*h;
	vector er =(k12-k0)*h;
	return (yh, er);
}

}

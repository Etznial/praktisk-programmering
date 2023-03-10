using System;
public static partial class ODE{

public static (vector,vector) rkstep23
        (Func<double,vector,vector> F, double x, vector y, double h)
{// Embedded Runge-Kutta stepper of the order 2/3
        vector k0 = F(x,y);  
        vector k1 = F(x+h/2  , y+(h/2  )*k0);
        vector k2 = F(x+3*h/4, y+(3*h/4)*k1);
        vector ka = (2*k0+3*k1+4*k2)/9;
        vector kb = k1;
        vector yh = y+ka*h;
        vector er = (ka-kb)*h;
        return (yh,er);
}//rkstep23

}//class

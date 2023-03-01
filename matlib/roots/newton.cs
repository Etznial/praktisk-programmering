using System;
using static System.Console;
public partial class roots{

public static vector newton
(Func<vector,vector> f, vector x, double eps=1e-3, vector dx=null){
	int n=x.size;
	vector fx=f(x),z,fz;
	while(true){
		matrix J=jacobian(f,x,fx,dx);
		var qrJ=new GivensQR(J);
		vector Dx=qrJ.solve(-fx);
		double s=1;
		while(true){// backtracking linesearch
			z=x+Dx*s;
			fz=f(z);
			if(fz.norm()<(1-s/2)*fx.norm()){ break; }
			if(s<1.0/64){ break; }
			s/=2;
		}
		x=z;
		fx=fz;
		if(fx.norm()<eps)break;
	}
	return x;
}//newton

}//roots

using System;
using static System.Console;
public partial class roots{

public static vector broyden
(Func<vector,vector> f, vector x, double eps=1e-3){
	vector fx=f(x),z,fz;
	matrix J = jacobian(f,x,fx);
	var qrj = new GivensQR(J);
	matrix B=qrj.inverse();
	while(true){
		vector Dx=-B*fx;
		double s=1;
		while(true){
			z=x+Dx*s;
			fz=f(z);
			if(fz.norm()<(1-s/2)*fx.norm()){
				break;
				}
			if(s<1.0/32){
				J=jacobian(f,z,fz);
				qrj=new GivensQR(J);
				B=qrj.inverse();
				break;
				}
			s/=2;
		}
		vector dx=z-x;
		vector df=fz-fx;

		vector c=(dx-B*df)/dx.dot(df); B.update(c,dx);
		//vector c=(dx-B*df)/(df%df); B.update(c,df);
		//vector c=(dx-B*df)/(dx%(B*df)); B.update(c,B*dx);
		
		x=z;
		fx=fz;
		if(fx.norm()<eps)break;
	}
	return x;
}//broyden

}//class

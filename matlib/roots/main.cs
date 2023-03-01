using System;
using static System.Math;
using static System.Console;
class main{
static void Main(){

vector c = new vector("1,1");
int n=c.size;
matrix A = new matrix(n,n);
var rnd=new System.Random(1);
for(int i=0;i<n;i++)
for(int j=0;j<n;j++)A[i,j]=2*(rnd.NextDouble()-0.5);
int ncalls=0;
Func<vector,vector> f;

f = delegate(vector z){
	ncalls++;
	vector r=new vector(2);
	double x=z[0],y=z[1],b=100;
	r[0]=2*(1-x)*(-1)+b*2*(y-x*x)*(-1)*2*x;
	r[1]=b*2*(y-x*x);
	return r;
	};

f = z => {// Rosenbrock
	ncalls++;
	vector r=new vector(2);
	double x=z[0],y=z[1],a=1,b=100;
	r[0]=a*(1-x);
	r[1]=b*(y-x*x);
	return r;
	};

//f = (z)=>{ ncalls++; return (A*(z-c)).map(t=>t*t*t); };
//f = (z)=>{ ncalls++; return (A*(z-c)).map(t=>Sin(t)); };
//double ene=3;
//f = (z)=>{ ncalls++; return new vector(ene*z[0]-z[1]*z[0],z[0]*z[0]-1); };

double eps=1e-7;
vector start = new vector("7 -7");

start.print("broyden:\nstart point:");
	ncalls=0; vector p = roots.broyden(f,start,eps);
	WriteLine($"ncalls={ncalls}");
	p.print("root=");
	f(p).print("f(root)=");
	WriteLine($"eps            = {eps}");
	WriteLine($"f(root).norm() = {f(p).norm()}");
	if(f(p).norm()<eps)WriteLine("test passed");
	else                  WriteLine("test failed");

start.print("\nnewton:\nstart point:");
	ncalls=0; p = roots.newton(f,start,eps);
	WriteLine($"ncalls={ncalls}");
	p.print("root=");
	f(p).print("f(root)=");
	WriteLine($"eps            = {eps}");
	WriteLine($"f(root).norm() = {f(p).norm()}");
	if(f(p).norm()<eps)WriteLine("test passed");
	else                  WriteLine("test failed");
}
}

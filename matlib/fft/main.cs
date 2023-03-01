using static System.Console;
using static complex;
using static cmath;
using static System.Math;
static class main{

static void Main(string[] argv){
if(argv.Length>0){
	int N=int.Parse(argv[0]);
	WriteLine($"N={N}");
	var xs = new complex[N];
	for(int i=0;i<N;i++) xs[i]=I;
	matlib.fft(xs);
	}
else	{
	var RND = new System.Random();
	//int N=128;
	int N=3*64;
	double noise=3,Afilter=0.2,Ffilter=0.025;
	var ts=new complex[N];
	var xs=new complex[N];
	for(int i=0;i<N;i++){
		ts[i]=(double)i/N;
		double omega=2*PI;
		xs[i]=2*sin(omega*ts[i])+cos(2*omega*ts[i]-2);
		}
	for(int i=0;i<N;i++) WriteLine($"{ts[i].Re} {xs[i].Re}");
	Write("\n\n");
	for(int i=0;i<N;i++)
		xs[i] += noise*RND.NextDouble()*(RND.NextDouble()-0.5)*2;
	for(int i=0;i<N;i++) WriteLine($"{ts[i].Re} {xs[i].Re}");
	Write("\n\n");
	var cs = matlib.fft(xs);
	double cmax=0;
	for(int i=0;i<N;i++) cmax = Max(cmax,abs(cs[i]));
	var cA=new complex[N];
	for(int i=0;i<N;i++)
		if(abs(cs[i]) < cmax*Afilter) cA[i]=0; else cA[i]=cs[i];
	var ys=matlib.ift(cA);
	for(int i=0;i<N;i++) WriteLine($"{ts[i].Re} {ys[i].Re}");
	Write("\n\n");
	var cF = new complex[N];
	for(int i=0;i<N;i++)
		if(i>N*Ffilter && i<N-N*Ffilter) cF[i]=0; else cF[i]=cs[i];
	ys = matlib.ift(cF);
	for(int i=0;i<N;i++) WriteLine($"{ts[i].Re} {ys[i].Re}");
	}
}//Main
}//class

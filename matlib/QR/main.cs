using static System.Console;
class main{
static void Main(string[] args){
int n=4,m=n;
if(args.Length>0) n=(int)double.Parse(args[0]);
if(args.Length>1) m=(int)double.Parse(args[1]);
var A=new matrix(n,m);
var rnd=new System.Random();
for(int i=0;i<A.size1;i++)
	for(int j=0;j<A.size2;j++)
		A[i,j]=2*(rnd.NextDouble()-1);
A.print("random matrix A:");
var qra=new GivensQR(A);
var qrb=new GSQR(A);
if(n>m){
	var pinvA=qrb.pinverse();
	pinvA.print("\npinvA (pseudo inverse of A):");
	(pinvA*A).print("\npinvA*A:");
	return;
}
var b=new vector(m);
for(int i=0;i<b.size;i++)
	b[i]=rnd.NextDouble();
b.print("random vector b:\n");
var x=qra.solve(b);
var y=qrb.solve(b);
x.print("GivnesQR: solution x to system A*x=b:\n");
y.print("GSQR    : solution x to system A*x=b:\n");
var Ax=A*x;
Ax.print("check: A*x (should be equal b):\n");
if(vector.approx(Ax,b)) WriteLine("test passed");
else WriteLine("test failed");
var B=qra.inverse();
var C=qrb.inverse();
(B*A).print("Givens:\nA^-1*A=");
(A*B).print("Givens:\nA*A^-1=");
(C*A).print("GramSchmidt:\nA^-1*A=");
(A*C).print("GramSchmidt:\nA*A^-1=");
}
}

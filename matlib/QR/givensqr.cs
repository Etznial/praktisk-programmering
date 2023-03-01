using static System.Math;
public class GivensQR{

public matrix QR;

public GivensQR(matrix M){
matrix A = M.copy();
for(int q=0;q<A.size2;q++){
	for(int p=q+1;p<A.size1;p++){
      		double theta=Atan2(A[p,q],A[q,q]);
		double c=Cos(theta),s=Sin(theta);
		for(int k=q;k<A.size2;k++){
			double xq=A[q,k], xp=A[p,k];
			A[q,k]= xq*c+xp*s;
			A[p,k]=-xq*s+xp*c;
			}
		A[p,q]=theta;
		}
	}
QR = A;
}

public vector solve(vector r){
	vector b=r.copy();
	for(int q=0;q<QR.size2;q++){
		for(int p=q+1;p<QR.size1;p++){
			double theta = QR[p,q];
			double c=Cos(theta),s=Sin(theta);
			double bq=b[q], bp=b[p];
			b[q]=+bq*c+bp*s;
			b[p]=-bq*s+bp*c;
			}
		}
	vector x = new vector(QR.size2);
	for(int i=QR.size2-1;i>=0;i--){
		double s=0; for(int k=i+1;k<QR.size2;k++) s+=QR[i,k]*x[k];
		x[i]=(b[i]-s)/QR[i,i];
		}
	return x;
	}//solve

public matrix inverse(){
	int m=QR.size2;
	var B=new matrix(m,m);
	var e=new vector(m);
	for(int i=0;i<m;i++){
		e[i]=1;
		B[i]=solve(e);
		e[i]=0;
		}
	return B;
	}//inverse
}

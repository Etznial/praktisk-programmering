using System;
public class GSQR{

public matrix Q,R;

public GSQR(matrix A){
	int m=A.size2;
	Q=A.copy();
	R=new matrix(m,m);
	for(int i=0;i<m;i++){
		R[i,i]=Q[i].norm();
		Q[i]/=R[i,i];
		for(int j=i+1;j<m;j++){
			R[i,j]=Q[i].dot(Q[j]);
			Q[j]-=Q[i]*R[i,j];
			}
		}
	}

public vector solve(vector b, bool useQ=true){
	vector x;
	if(useQ) x=Q%b;
	else x = b.copy();
	for(int i=x.size-1;i>=0;i--){
		double s=0;for(int k=i+1;k<x.size;k++)s+=R[i,k]*x[k];
		x[i]=(x[i]-s)/R[i,i];
		}
	return x;
	}

public matrix pinverse(){
        int m=R.size2;
        var B=new matrix(m,m);
        var e=new vector(m);
        for(int i=0;i<m;i++){
                e[i]=1;
                B[i]=solve(e,useQ:false);
                e[i]=0;
                }
        return B*Q.T;
        }//pinverse

public matrix inverse(){
	if(Q.size1!=Q.size2)throw new Exception("GSQR: inverse: size1!=size2");
        int m=R.size2;
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

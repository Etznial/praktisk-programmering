using static System.Console;
using static System.Math;
public class vec{
	public double x,y,z;
	public vec (double a, double b, double c){x=a;y=b;z=c;}
	public void print(string s){Write(s);WriteLine($"{x} {y} {z}");}
	public void print(){this.print("");}
	
	public override string ToString(){return $"({x} {y} {z})";}

	public static vec operator+(vec u, vec v){/* u+v */
		return new vec(u.x+v.x, u.y+v.y, u.z+v.z);
			} 
	public static vec operator-(vec u, vec v){/* u-v */
		return new vec(u.x-v.x, u.y-v.y, u.z-u.z);
	}

	public static vec operator-(vec u){/* -u */
		return new vec(-u.x, -u.y, -u.z);
	}	


	public static vec operator*(vec u, double c){/* u*c */
		return new vec(u.x*c,u.y*c,u.z*c);
	}
	
	public static vec operator*(double c, vec u){ /* c,*u */
		return u*c;
	}
	
	public static double operator% (vec u, vec v){ /* u%v => dot product */
		return u.x*v.x + u.y*v.x + u.z*v.z;
	}

	public double dot (vec other){ return this%other;}

	public vec cross(vec other){
		return new vec(this.y*other.z-other.y*this.z, this.z*other.x-other.z*this.x, this.x*other.y-other.x*this.y);
	} 

	public double norm(){
		return Sqrt(x*x+y*y+z*z);
	}

}

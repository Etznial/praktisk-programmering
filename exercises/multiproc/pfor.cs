using System;
using System.Threading;
using System.Threading.Tasks;
using static System.Console;
using static System.Math;

class main{

public static int Main(string[] args){
	int nterms=(int)1e8, nthreads=1;	
	foreach(var arg in args)
	{
		var words = arg.Split(':');
		if(words[0]=="-terms"    )nterms=(int)float.Parse(words[1]);
		if(words[0]=="-threads")nthreads=(int)float.Parse(words[1]);
	}	
	
	WriteLine("=======================[new threading with Parallel.For]=======================");
	double sum = 0;
	//for(int i=1;i<nterms+1;i++){sum+=1.0/i;}
	Parallel.For( 1, nterms+1, delegate(int i) { sum+=1.0/i;} );
	WriteLine($"nterms = {nterms} nthreads = {nthreads} sum = {sum}");
	WriteLine($"Grunden til at det er hurtigerer at køre et serial for-loop end det er for Parallel.For er fordi man kalder delegate(int i) i Parallel.For mega mange gange, hvor man ikke kalder en function i det simple for-loop, har man meget mere komplekse funktioner end at summere tal, så kan man godt gå tilbage til Parallel.For fordi så er tiden det tager at kalde delegate() meget lille i forhold til den mere komplekse funktion man multithreader. Vil man undgå at kalde delegate skal man multithreade manuelt");
	return 0;
	}
}

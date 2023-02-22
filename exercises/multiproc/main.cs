using System;
using System.Threading;
using static System.Console;
using static System.Math;
class main{

public class data {public int a,b; public double sum;}

public static void harmonic(object obj){
	var local = (data)obj;
	local.sum=0;
	for(int i=local.a;i<local.b;i++)local.sum+=1.0/i;
	WriteLine($"{Thread.CurrentThread.Name} a={local.a} b={local.b}");
	WriteLine($"{Thread.CurrentThread.Name} partial sum = {local.sum}");
	}


public static int Main(string[] args){
	int nterms=(int)1e8, nthreads=1; /*Start værdier*/	
	foreach(var arg in args)
	{
		var words = arg.Split(':');
		if(words[0]=="-terms"  )nterms  =(int)float.Parse(words[1]);
		if(words[0]=="-threads")nthreads=(int)float.Parse(words[1]);
	}
	WriteLine("=======================[new threading]=======================");
	WriteLine($"nterms = {nterms} nthreads = {nthreads}");

	data[] x = new data[nthreads]; /*Separate threds*/
	for(int i=0;i<nthreads;i++){
		x[i] = new data();
		x[i].a = 1 + i*(nterms/nthreads); /*Paranteser skal flyttes for ikke at få integer overflow*/
		x[i].b = 1 + (i+1)*(nterms/nthreads); /*Paranteser skal flyttes for ikke at få integer overflow*/
		WriteLine($"x.a = {x[i].a} x.b = {x[i].b}");
		}
	x[x.Length-1].b=nterms+1; /*Check endpoint*/

	Thread[] threads = new Thread[nthreads];
	for(int i=0;i<nthreads;i++){
	 	threads[i] = new Thread(harmonic); /*Create a new thread*/
		threads[i].Name = $"thread number {i+1}";  /* Changes the name of the threads for better user reading*/
		threads[i].Start(x[i]);
		}
	
	WriteLine("master thread: now waiting for other threads...");
	for(int i=0;i<nthreads;i++){ /* Joining the threads */
		threads[i].Join();
	}
	double total = 0;
	for(int i=0;i<nthreads;i++){
		total+=x[i].sum;
	}
	WriteLine($"total sum = {total}");	
	return 0;
	}
}

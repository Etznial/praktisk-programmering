using System;
using static System.Console;
class main{
public static void Main(string[] argv){
	int N=100;
	if(argv.Length>0)N=int.Parse(argv[0]);
	var qnd = new halton(dim:2,shift:1);
	var rnd = new Random();
	for(int i=0;i<N;i++){
		var x=qnd.next();
		WriteLine($"{x[0]} {x[1]}");
		Error.WriteLine($"{rnd.NextDouble()} {rnd.NextDouble()}");
		}
}
}//main

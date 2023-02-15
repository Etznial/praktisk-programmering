using System;
using static System.Console;
using static System.Math;
class main{
	public static int Main(string[] args){
	WriteLine($"Input: {args}");
	foreach(var arg in args){
		var words = arg.Split(':');
		if(words[0]=="-numbers"){
			var numbers=words[1].Split(',');
			foreach(var number in numbers){
				double x = double.Parse(number);
				WriteLine($"{x} {Sin(x)} {Cos(x)}");
				}
			}
		}
	char[] split_delimiters = {' ','\t','\n'};
	var split_options =StringSplitOptions.RemoveEmptyEntries;
	for( string line = ReadLine(); line != null; line = ReadLine()){
		WriteLine($"Input: {line}");
		var numbers = line.Split(split_delimiters,split_options);
		foreach(var number in numbers){
			double x = double.Parse(number);
				WriteLine($"{x} {Sin(x)} {Cos(x)}");	
			}
		}
	string infile=null,outfile=null;
	foreach(var arg in args){
		var words=arg.Split(':');
		if(words[0]=="-input")infile=words[1];
		if(words[0]=="-output")outfile=words[1];
		}
	if(infile==null || outfile==null){
		Error.WriteLine("wrong filename argument");
		return 1;
		}
	var instream =new System.IO.StreamReader(infile);
	var outstream=new System.IO.StreamWriter(outfile,append:false);
	for(string line=instream.ReadLine();line!=null;line=instream.ReadLine()){
		double x = double.Parse(line);
		outstream.WriteLine($"{x} {Sin(x)} {Cos(x)}");
		}
	instream.Close();
	outstream.Close();
	return 0;
	}
}

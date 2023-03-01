using static System.Console;

class main{
static void Main(string[] args){
	foreach(var arg in args){
		var words = arg.Split(':');
		if(words[0] == "-type"){
			if(words[1] == "gamma"){
				for(double x=-5+1.0/128;x<=5;x+=1.0/64){
					WriteLine($"{x} {sfuns.gamma(x)}");
					}
				}
			if(words[1] == "error"){
		 		for(double x=-4+1.0/128;x<=4;x+=1.0/64){
					WriteLine($"{x} {sfuns.erf(x)}");
					}
				} 
		}
	}
}

}// class

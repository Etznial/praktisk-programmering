Out.txt : main.exe
	mono main.exe > Out.txt

main.exe : main.cs matrix.dll
		mcs main.cs -reference:matrix.dll -target:exe -out:main.exe

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/matrix.cs $(DIR)/vector.cs
	mcs $< -target:library -out:$@

clean : 
	rm -f main.exe matrix.dll Out.txt

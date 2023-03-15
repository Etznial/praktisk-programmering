CSC = mcs


Out.txt : main.exe
	mono main.exe -rmax:10 -N:11 > Out.txt

main.exe : main.cs matrix.dll jacobi.dll
	$(CSC) main.cs -reference:matrix.dll -reference:jacobi.dll -target:exe -out:main.exe

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	$(CSC) -target:library -out:$@ $^

jacobi.dll : jacobi.cs matrix.dll
	$(CSC) -target:library -reference:matrix.dll -out:$@ jacobi.cs

clean : 
	$(RM) *.dll *.exe [Oo]ut*

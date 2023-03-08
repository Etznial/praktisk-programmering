CSC = mcs


Out.txt : main.exe
	mono main.exe > Out.txt

main.exe : main.cs matrix.dll
	$(CSC) main.cs -reference:matrix.dll -target:exe -out:main.exe

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	$(CSC) -target:library -out:$@ $^

clean : 
	$(RM) *.dll *.exe [Oo]ut*

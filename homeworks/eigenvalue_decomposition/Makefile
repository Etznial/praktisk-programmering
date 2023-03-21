default : main.exe Out.txt Eigenrmax.svg Eigendr.svg Eigenfunc.svg

CSC = mcs

Out.txt : main.exe 
	mono main.exe -rmax:9 -dr:0.1 > Out.txt

main.exe : main.cs matrix.dll jacobi.dll
	$(CSC) main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

Eigenrmax.svg : eigenrmax.gpi main.exe rmax_lowest.data
	gnuplot $<


Eigendr.svg : eigendr.gpi main.exe dr_lowest.data
	gnuplot $<

Eigenfunc.svg : eigenfunc.gpi main.exe eigenfunc0.data eigenfunc1.data eigenfunc2.data
	gnuplot $<


DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	$(CSC) -target:library -out:$@ $^

jacobi.dll : jacobi.cs matrix.dll
	$(CSC) -target:library -reference:matrix.dll -out:$@ jacobi.cs

clean : 
	$(RM) *.data *.dll *.exe [Oo]ut*

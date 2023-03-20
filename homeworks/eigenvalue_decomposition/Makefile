CSC = mcs

Out.txt : main.exe Eigenrmax.svg Eigendr.svg Eigenfunc.svg
	mono main.exe -rmax:9 -dr:0.1 > Out.txt

Eigenrmax.svg : eigenrmax.gpi main.exe rmax_lowest.data
	gnuplot $<

Eigendr.svg : eigendr.gpi main.exe dr_lowest.data
	gnuplot $<

Eigenfunc.svg : eigenfunc.gpi main.exe eigenfunc.data
	gnuplot $<

main.exe : main.cs matrix.dll jacobi.dll
	$(CSC) main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	$(CSC) -target:library -out:$@ $^

jacobi.dll : jacobi.cs matrix.dll
	$(CSC) -target:library -reference:matrix.dll -out:$@ jacobi.cs

clean : 
	$(RM) *.dll *.exe [Oo]ut*

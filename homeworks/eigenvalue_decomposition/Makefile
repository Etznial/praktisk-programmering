default : Eigenrmax.svg Eigendr.svg Eigenfunc.svg

Out.txt eigenfunc0.data eigenfunc1.data eigenfunc2.data rmax_lowest.data dr_lowest.data: main.exe 
	mono main.exe -rmax:9 -dr:0.1 > Out.txt -target:$@

main.exe : main.cs matrix.dll jacobi.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

Eigenrmax.svg : eigenrmax.gpi Out.txt rmax_lowest.data
	gnuplot $<


Eigendr.svg : eigendr.gpi Out.txt dr_lowest.data
	gnuplot $<

Eigenfunc.svg : eigenfunc.gpi Out.txt eigenfunc0.data eigenfunc1.data eigenfunc2.data
	gnuplot $<


DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	mcs -target:library -out:$@ $^

jacobi.dll : jacobi.cs matrix.dll
	mcs -target:library -reference:matrix.dll -out:$@ jacobi.cs

clean : 
	$(RM) *.data *.dll *.exe [Oo]ut*

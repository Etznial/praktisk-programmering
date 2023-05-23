default : Out.txt

Out.txt : main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll genlist.dll QRGS.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

DIR1 = ../../matlib/matrix
matrix.dll : $(DIR1)/vector.cs $(DIR1)/matrix.cs
	mcs -target:library -out:$@ $^

DIR2 = ../../exercises/genlist/
genlist.dll : $(DIR2)/genlist.cs
	mcs -target:library -out:$@ $^

DIR3 = ../linear_equations/
QRGS.dll : $(DIR3)QRGS.cs matrix.dll
	mcs $(addprefix -reference:,$(filter %.dll,$^)) -target:library -out:$@ $(DIR3)QRGS.cs

PenBtest.svg : penBtest.gpi penB.data pen.data
	gnuplot $<


#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@

clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data *.svg

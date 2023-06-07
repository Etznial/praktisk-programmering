default : Out.txt

Out.txt : main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll genlist.dll sfuns.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

run : 
	python3.10 mainpy.py 

DIR1 = ../../matlib/matrix
matrix.dll : $(DIR1)/vector.cs $(DIR1)/matrix.cs
	mcs -target:library -out:$@ $^

DIR2 = ../../exercises/genlist/
genlist.dll : $(DIR2)/genlist.cs
	mcs -target:library -out:$@ $^

DIR3 = ../../exercises/plots/
sfuns.dll : $(DIR3)error.cs $(DIR3)gamma.cs
	mcs -target:library -out:$@ $^

Erf.svg : erf.gpi erf.data erfTable.data
	gnuplot $<

erfTable.data : Makefile
	>$@
	echo 0 0 >>$@
	echo 1 0.842700793 >>$@
	echo 2 0.995322265 >>$@
	echo 3 0.999977910 >>$@

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data *.svg

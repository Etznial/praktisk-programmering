default : Out.txt Higgs.svg

Out.txt bwfit.data : main.exe 
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll genlist.dll sfuns.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

DIR1 = ../../matlib/matrix
matrix.dll : $(DIR1)/vector.cs $(DIR1)/matrix.cs
	mcs -target:library -out:$@ $^

DIR2 = ../../exercises/genlist/
genlist.dll : $(DIR2)/genlist.cs
	mcs -target:library -out:$@ $^

DIR3 = ../../exercises/plots/
sfuns.dll : $(DIR3)error.cs $(DIR3)gamma.cs
	mcs -target:library -out:$@ $^

Higgs.svg : higgs.gpi higgs.data bwfit.data
	gnuplot $<

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.svg bwfit.data

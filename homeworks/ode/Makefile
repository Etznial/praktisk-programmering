default : Out.txt Rktest.svg Drivertest.svg Pentest.svg

Out.txt rk.data driver.data: main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll genlist.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

DIR1 = ../../matlib/matrix
matrix.dll : $(DIR1)/vector.cs $(DIR1)/matrix.cs
	mcs -target:library -out:$@ $^

DIR2 = ../../exercises/genlist/
genlist.dll : $(DIR2)/genlist.cs
	mcs -target:library -out:$@ $^
	

Rktest.svg : rktest.gpi rk.data
	gnuplot $<

Drivertest.svg : drivertest.gpi driver.data
	gnuplot $<

Pentest.svg : pentest.gpi pen.data
	gnuplot $<

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

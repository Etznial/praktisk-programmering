default : Out.txt Rktest.svg Drivertest.svg Pentest.svg PenBtest.svg Utest.svg

Out.txt rk.data driver.data pen.data pen1.data penB.data u1.data u2.data u3.data : main.exe mainB.exe
	mono main.exe > Out.txt -target:$@
	mono mainB.exe >> Out.txt -target:$@

main.exe : main.cs matrix.dll genlist.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

mainB.exe : mainB.cs matrix.dll genlist.dll
	mcs mainB.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:mainB.exe

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

Pentest.svg : pentest.gpi pen.data pen1.data
	gnuplot $<

PenBtest.svg : penBtest.gpi penB.data pen.data
	gnuplot $<

Utest.svg : utest.gpi u1.data u2.data u3.data
	gnuplot $<

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

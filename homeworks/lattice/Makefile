default : Out.txt Lattice.svg

Out.txt testLattice1.data testLattice2.data testLattice3.data : main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

Lattice.svg : lattice.gpi testLattice1.data testLattice2.data testLattice3.data
	gnuplot $<


DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	mcs -target:library -out:$@ $^

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

default : Out.txt Akima.svg

Out.txt akima.data: main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

Akima.svg : akima.gpi akima.data
	gnuplot $<


DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	mcs -target:library -out:$@ $^

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

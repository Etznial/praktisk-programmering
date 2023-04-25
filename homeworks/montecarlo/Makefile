default : Out.txt Testplainmc.svg Testqrmc.svg

Out.txt coscosIt.data QRMCcoscosIt.data : main.exe
	mono main.exe > Out.txt -target:$@

main.exe : main.cs matrix.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	mcs -target:library -out:$@ $^

Testplainmc.svg : testplainmc.gpi coscosIt.data
	gnuplot $<

Testqrmc.svg : testqrmc.gpi QRMCcoscosIt.data coscosIt.data
	gnuplot $<

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

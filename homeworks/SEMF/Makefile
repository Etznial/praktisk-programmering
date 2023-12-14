default : Out.txt BindingEnergy.svg MassEnergy.svg MassEnergypy.svg

Out.txt BE37Cl.data ME37Cl.data ME37Clpy.data : main.exe semf.py
	mono main.exe > Out.txt -target:$@
	python3 semf.py >> Out.txt -target:$@

main.exe : main.cs matrix.dll
	mcs main.cs $(addprefix -reference:,$(filter %.dll,$^)) -target:exe -out:main.exe

BindingEnergy.svg : bindingEnergy.gpi BE37Cl.data
	gnuplot $<

MassEnergy.svg : massEnergy.gpi ME37Cl.data
	gnuplot $<

MassEnergypy.svg : massEnergypy.gpi ME37Clpy.data
	gnuplot $<

DIR = ../../matlib/matrix
matrix.dll : $(DIR)/vector.cs $(DIR)/matrix.cs
	mcs -target:library -out:$@ $^

#QRGS.dll : QRGS.cs matrix.dll
#	mcs QRGS.cs -target:library -reference:matrix.dll -out:$@
clean : 
	$(RM) *.dll *.exe [Oo]ut* *.data

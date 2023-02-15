all: Out.txt

Out.txt : main.exe
	mono main.exe > Out.txt  2> log 

main.exe : main.cs cmath.dll 
	mcs -target:exe -out:$@ \
		$(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

cmath.dll : cmath.cs complex.cs
	mcs -target:library -out:cmath.dll $^

.PHONEY : clean
clean:
	$(RM) *.exe *.dll [Oo]ut* log*

all: Out.txt
	cat $<
Out.txt: main.exe
	mono $< > $@
main.exe: main.cs cmath.dll
	mcs $(filter %.cs,$^) -target:exe -reference:$(filter %.dll,$^) -out:$@
cmath.dll: cmath.cs complex.cs 
	mcs $^ -target:library -out:$@
clean:
	$(RM) [Oo]ut* *.dll *.exe

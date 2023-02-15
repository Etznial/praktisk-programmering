#RM = rm ---force

all: Out.txt

Out.txt : main.exe Makefile infile.txt numbers.txt
	mono main.exe hello world \
		-input:infile.txt \
		-output:outfile.txt \
		-numbers:1,2,3,4,5,6> $@ 2> log < infile.txt
	cat infile.txt | mono main.exe hello world \
		-input:infile.txt \
		-output:outfile.txt \
		-numbers:1,2,3,4,5,6 1>> $@ 2>> log
	cat infile.txt | mono main.exe hello world \
		-input:infile.txt \
		-output:outfile.txt \
		$$(cat numbers.txt) 1>> $@ 2>> log

main.exe : main.cs input.dll
	mcs -target:exe -out:$@ \
		$(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

input.dll : input.cs
	mcs -target:library -out:$@ $(filter %.cs,$^)

infile.txt : Makefile
	echo 1 >>$@
	echo 2 >>$@
	echo 3 >>$@

numbers.txt:Makefile
	echo "-numbers:9,8,7" > $@

.PHONEY : clean
clean: 
	$(RM) *.exe *.dll [Oo]ut* log*

test: 
	echo $(RM)

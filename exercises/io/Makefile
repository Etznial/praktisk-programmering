all: Out.txt

Out.txt : main.exe Makefile numbers.txt input.txt input2.txt
	mono main.exe -numbers:1,2,3,4,5 > Out.txt  2> log \
		-output:output.txt < input.txt > $@ \
		-input:input2.txt \

main.exe : main.cs
	mcs -target:exe -out:$@ \
		$(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

input.txt : Makefile
	echo 6 7 8 9 10 > input.txt

input2.txt : Makefile
	echo "11\n12\n13\n14\n15" > input2.txt

.PHONEY : clean
clean:
	$(RM) *.exe *.dll [Oo]ut* log*

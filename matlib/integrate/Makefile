all:Out.htm
	cat $<

Out.htm: main.exe test.py Makefile
	echo '<meta charset="UTF-8"><pre>' >$@
	mono $< >> $@ 2>log
	python3 test.py >> $@

main.exe:main.cs math.dll
	mcs $< -reference:math.dll

math.dll: integrate.cs adapt3.cs
	mcs $^ -target:library -out:$@

clean:
	rm -f *.exe *.dll [Oo]ut* log*

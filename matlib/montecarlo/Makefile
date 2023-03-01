PLOT = Points.png
$(PLOT): out.quasi.txt out.pseudo.txt Makefile
	echo '\
	set term $(subst .,,$(suffix $@));\
	set out "$@";\
	set tics out;\
	set xtics 0.2;\
	set ytics 0.2;\
	set multiplot;\
	set size 0.5,0.5;\
	set origin 0.0,0.5;\
	plot "out.pseudo.txt" with points pointtype 5 pointsize 0.2 notitle;\
	set size 0.5,0.5;\
	set origin 0.5,0.5;\
	plot "out.quasi.txt" with points pointtype 5 pointsize 0.2 notitle;\
	'|gnuplot

N = 1000
out.quasi.txt out.pseudo.txt: main.exe Makefile
	mono $< $(N) 1> out.quasi.txt 2>out.pseudo.txt

main.exe: main.cs matlib.dll
	mcs -t:exe -out:$@ -r:matlib.dll main.cs

LIBSRC = ../matrix/vector.cs
matlib.dll: $(LIBSRC) halton.cs
	mcs -target:library -out:$@ $^

clean:
	$(RM) [Oo]ut* *.exe *.dll *.png *.svg *.gif

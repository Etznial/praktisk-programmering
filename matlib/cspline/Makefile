FMT = svg
OPT = font "Arial,20" size 1024,600
PLOT = Cspline
PLOTFILE=$(PLOT).$(FMT)
DATA = out.txt
ECHO=/bin/echo

$(PLOTFILE): $(PLOT).gpi $(DATA)
	gnuplot $<

$(DATA): main.exe
	mono $< >$@

main.exe: main.cs matlib.dll
	mcs $< -reference:$(lastword $^)

matlib.dll: cspline.cs
	mcs $< -target:library -out:$@

$(PLOT).gpi: Makefile
	cat /dev/null >$@
	$(ECHO) 'set term $(FMT) $(OPT)' >>$@
	$(ECHO) 'set out "$(PLOTFILE)"' >>$@
	$(ECHO) 'set multiplot layout 2,2 title "cspline, derivative, integral"' >>$@
	$(ECHO) 'set key top right' >>$@
	$(ECHO) 'set xtics 1' >>$@
	$(ECHO) 'set grid' >>$@
	$(ECHO) 'set rmargin 0' >>$@
	$(ECHO) 'plot [][-1:1]\' >>$@
	$(ECHO) ' "$(DATA)" index 0 with points pointtype 7 title "points"\' >>$@
	$(ECHO) ',"$(DATA)" index 1 u 1:2 with lines title "exact"\' >>$@
	$(ECHO) ',"$(DATA)" index 1 u 1:3 with lines title "cspline"' >>$@
	$(ECHO) 'set key center top' >>$@
	$(ECHO) 'set ytics format ""' >>$@
	$(ECHO) 'set rmargin' >>$@
	$(ECHO) 'plot [][-1:1]\' >>$@
	$(ECHO) ' "$(DATA)" index 2 u 1:2 with lines title "exact derivative"\' >>$@
	$(ECHO) ',"$(DATA)" index 2 u 1:3 with lines title "spline derivative"' >>$@
	$(ECHO) 'set ytics format "%h"' >>$@
	$(ECHO) 'set key center bottom' >>$@
	$(ECHO) 'set rmargin 0' >>$@
	$(ECHO) 'plot [][0:2]\' >>$@
	$(ECHO) ' "$(DATA)" index 3 u 1:2 with lines title "exact integral"\' >>$@
	$(ECHO) ',"$(DATA)" index 3 u 1:3 with lines title "spline integral"' >>$@
	$(ECHO) 'unset multiplot' >>$@

clean:
	rm -f *.exe *.dll $(PLOTFILE) $(DATA) *.gpi *.svg out*

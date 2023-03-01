SDKDIR=$(HOME)/.dotnet/sdk/6.0.102# or whatever your sdk is
CSCPATH=$(SDKDIR)/Roslyn/bincore/csc.dll
NETSTANDARD=$(SDKDIR)/ref/netstandard.dll
CONFIG=$(SDKDIR)/vstest.console.runtimeconfig.json
DOTNET=DOTNET_CLI_TELEMETRY_OPTOUT=1; dotnet
CSC=$(DOTNET) $(CSCPATH) -reference:$(NETSTANDARD)
RUN=$(DOTNET) exec --runtimeconfig $(CONFIG)
CSC = mcs
RUN = mono

DATA-SIN = out.sin.data.txt
DATA-ERF = out.erf.data.txt
PLOT-SIN = Plot.sin.png
PLOT-ERF = Plot.erf.png
all: $(PLOT-SIN) $(PLOT-ERF)

$(DATA-SIN): main.exe
	$(RUN) $< 1> $@ 2> log

$(DATA-ERF): main-erf.exe
	$(RUN) $< 1> $@ 2> log-erf

main.exe: main.cs matlib.dll
	$(CSC) $< -reference:matlib.dll

main-erf.exe: main-erf.cs erf.cs matlib.dll
	$(CSC) -out:$@ -reference:matlib.dll $(filter %.cs,$^)

matlib.dll: rkstep12.cs rkstep23.cs driver.cs ../matrix/vector.cs
	$(CSC) -target:library -out:$@ $^

clean:
	rm -f *.dll *.png *.exe out* *.svg [Ll]og*

$(PLOT-SIN): $(DATA-SIN) Makefile
	echo '\
	set term png;\
	set out "$(PLOT-SIN)";\
	set key out;\
	set tics out;\
	set grid;\
	set xlabel "x";\
	plot \
	 "$(DATA-SIN)" using 1:2 with linespoints pt 7 ps 0.5 title "sin" \
	,"$(DATA-SIN)" using 1:3 with linespoints pt 5 ps 0.5 title "cos" \
	,"$(DATA-SIN)" using ($$1):(0) w p pt 7 ps 0.1 not \
	'|gnuplot
	
$(PLOT-ERF): $(DATA-ERF) Makefile
	echo '\
	set term png;\
	set out "$@";\
	set key out;\
	set tics out;\
	set grid;\
	set xlabel "x";\
	set title "erf(x) as solution to differential equation";\
	plot \
	 "$<" index 0 using 1:2 with line title "erf(x)" \
	,"$<" index 1 using 1:2 with points pt 5 title "table" \
	'|gnuplot
	

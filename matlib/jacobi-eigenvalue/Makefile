SDKDIR=$(HOME)/.dotnet/sdk/6.0.102# or whatever your sdk is
CSCPATH=$(SDKDIR)/Roslyn/bincore/csc.dll
NETSTANDARD=$(SDKDIR)/ref/netstandard.dll
CONFIG=$(SDKDIR)/vstest.console.runtimeconfig.json
DOTNET=dotnet
CSC=$(DOTNET) $(CSCPATH) -reference:$(NETSTANDARD)
RUN=$(DOTNET) exec --runtimeconfig $(CONFIG)
# remember to export DOTNET_CLI_TELEMETRY_OPTOUT=1
CSC = mcs -optimize+ -platform:arm
RUN = mono -O=all #-O=all,-shared #--optimize=unsafe,loop,inline #--gc=sgen --llvm --optimize=all

CFLAGS = -Ofast -march=native
LDLIBS = $(shell gsl-config --libs)

comma:=,
empty:=
space:=$(empty) $(empty)
commalist = $(subst $(space),$(comma),$(1))

all: out.txt outJ.txt Times.svg Times.gsl.svg

out.txt: main.exe #Makefile
	$(RUN) $< 7 | tee $@

outJ.txt: mainJ.exe #Makefile
	$(RUN) $< 7 | tee $@

out.gsl.txt: main #Makefile
	cat /dev/null > $@
	for N in `seq 200 10 300`; do \
	echo "N=$$N";\time -ao $@ -f "$$N %U" ./main $$N ; \
	done

main:main.c

Times.gsl.svg: out.gsl.txt Makefile
	echo '\
	set term svg size 640,480 background "white" font "times,14";\
	set out "$@";\
	set key left;\
	set title "matrix diagonalization times (gsl `gsl-config --version` `uname -m`)";\
	set xlabel "matrix size n";\
	set ylabel "diagonalization time t, sec";\
	f(x)=b+(x/a)**c;\
	b=1; a=90; c=3;\
	fit f(x) "$<" via a,b;\
	plot \
	 "$<" title "measurement" \
	,f(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",a,c,b) \
	'| gnuplot

Times.svg: out.times.txt out.timesJ.txt out.times-nov.txt \
	out.times1.txt Makefile
	echo '\
	set term svg size 640,480 background "white" font "times,14";\
	set out "$@";\
	set key left;\
	set title "matrix diagonalization times (`$(CSC) --version` `uname -m`)";\
	set xlabel "matrix size n";\
	set ylabel "diagonalization time t, sec";\
	b=1; a=90; c=3;\
	f(x)=b+(x/a)**c;\
	bj=1; aj=90; cj=3;\
	fj(x)=bj+(x/aj)**cj;\
	bv=1; av=90; cv=3;\
	fv(x)=bv+(x/av)**cv;\
	b1=1; a1=90; c1=3;\
	f1(x)=b1+(x/a1)**c1;\
	fit f(x) "$<" via a,b;\
	fit fj(x) "out.timesJ.txt" via aj,bj;\
	fit fv(x) "out.times-nov.txt" via av,bv;\
	fit f1(x) "out.times1.txt" via a1,b1;\
	plot \
	 "out.timesJ.txt" title "timesJ/Jtimes with V" \
	,fj(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",aj,cj,bj) \
	,"$<" title "Jacobi with V" \
	, f(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",a,c,b) \
	,"out.times-nov.txt" title "Jacobi no V" \
	,fv(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",av,cv,bv) \
	,"out.times1.txt" title "Jacobi nvals=1" \
	,f1(x) title sprintf("fit: (n/%.f)^{%.2f}+%.3f",a1,c1,b1) \
	'| gnuplot

out.times.txt: main.exe
	cat /dev/null > $@
	for N in `seq 63 5 123`; do \
	echo "N=$$N";\time -ao $@ -f "$$N %U" $(RUN) $< $$N ; \
	done

out.times1.txt: main.exe
	cat /dev/null > $@
	for N in `seq 93 5 143`; do \
	echo "N=$$N";\time -ao $@ -f "$$N %U" $(RUN) $< $$N 1 ; \
	done

out.times-nov.txt: main-nov.exe
	cat /dev/null > $@
	for N in `seq 63 5 123`; do \
	echo "N=$$N";\time -ao $@ -f "$$N %U" $(RUN) $< $$N ; \
	done

out.timesJ.txt: mainJ.exe
	cat /dev/null > $@
	for N in `seq 53 5 103`; do \
	echo "N=$$N";\time -ao $@ -f "$$N %U" $(RUN) $< $$N ; \
	done

main.exe: main.cs matrix.dll
	$(CSC) $< -out:$@ -r:$(call commalist,$(filter-out $<,$^))

main-nov.exe: main-nov.cs matrix.dll
	$(CSC) $< -out:$@ -r:$(call commalist,$(filter-out $<,$^))

mainJ.exe: mainJ.cs matrix.dll
	$(CSC) $< -out:$@ $(addprefix -r:,$(filter %.dll,$^))

matrix.dll: EVD.cs jacobi.cs ../matrix/matrix.cs ../matrix/vector.cs
	$(CSC) -t:library -out:./$@ $^

clean: ; $(RM) *.dll *.exe out* *.log *.svg

#%.dll: %.cs ; $(CSC) -t:library $*.cs -o:$*.dll *log

test:
	echo $(LDLIBS)
	echo $$($(RUN) --version)

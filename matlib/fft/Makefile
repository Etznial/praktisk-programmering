all:Out.noise.svg Out.timing.svg

Out.noise.svg: out.noise.data Makefile
	echo '\
	set term svg size 640,480 background "white" font "Arial,18";\
	set out "$@";\
	set title "Noise filtering using DFT";\
	set xlabel "time";\
	set ylabel "signal";\
	plot \
	 "$<" index 0 with lines title "original signal"\
	,"$<" index 1 with lines title "signal with random noise"\
	,"$<" index 2 with lines title "amplitude filtered"\
	,"$<" index 3 with lines title "frequency filtered"\
	'|gnuplot

Out.timing.svg: out.timing.data Makefile
	echo '\
	set term svg size 640,480 background "white" font "Arial,18";\
	set out "$@";\
	set title "FFT time as function of data set size";\
	set key left;\
	set xlabel "data set size N";\
	set ylabel "seconds to perform FFT";\
	set format x "%3.0e";\
	a=1;\
	f(x)=a*x*log(x);\
	fit f(x) "$<"  via a;\
	plot "$<" with points title "data", f(x) with lines title "N*log(N)";\
	'|gnuplot

#N = 16384 32768 65536 131072 #262144 segmentation fault
#N = $(shell echo '2^14; 2^15; 2^16; 2^17'|bc -l)
N = $(shell echo '4*2^12; 4*2^13; 4*2^14; 4*2^15'|bc -l)
out.timing.data:main.exe #Makefile
	>$@
	for n in $(N); do time -f "$$n %U" -ao $@ mono ./$< $$n; done

out.noise.data:main.exe
	mono ./$< > $@

main.exe: main.cs fft.cs matlib.dll
	mcs -out:$@ $(addprefix -r:,$(filter %.dll,$^))  $(filter %.cs,$^)

matlib.dll: \
	../complex/complex.cs\
	../complex/cmath.cs
	mcs -target:library -out:$@ $^

clean:
	find . -type f -executable -delete
	$(RM) *.o [Oo]ut* *.log

TEST=$(shell echo '2^11; 2^12; 2^13'|bc -l)

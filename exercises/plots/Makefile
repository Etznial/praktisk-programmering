default : Gamma.svg Error.svg Lngamma.svg main.exe

Gamma.svg : gamma.gpi gamma.data factorials.data
	gnuplot $<

gamma.data : main.exe
	mono $< -type:gamma > $@

Error.svg : error.gpi error.data errorTable.data
	gnuplot $<

error.data : main.exe
	mono $< -type:error > $@

Lngamma.svg : lngamma.gpi lngamma.data lnfactorials.data
	gnuplot $<

lngamma.data : main.exe
	mono $< -type:lngamma > $@

main.exe : main.cs sfuns.dll
	mcs -target:exe -out:$@ \
		$(filter %.cs,$^) $(addprefix -reference:,$(filter %.dll,$^))

sfuns.dll : gamma.cs error.cs
	mcs -target:library -out:$@ $^

factorials.data : Makefile
	>$@
	echo 0 1 >>$@
	echo 1 1 >>$@
	echo 2 2 >>$@
	echo 3 6 >>$@
	echo 4 24 >>$@
	echo 5 120 >>$@
#	echo 6 $(echo "6*120" | bc) >>$@
#
lnfactorials.data : Makefile
	>$@
	echo 0 0 >>$@
	echo 1 0 >>$@
	echo 2 0.69315 >>$@
	echo 3 1.79175 >>$@
	echo 4 3.17805 >>$@
	echo 5 4.78749 >>$@


errorTable.data : Makefile
	>$@
	echo 0 0 >>$@
	echo 1 0.842700793 >>$@
	echo 2 0.995322265 >>$@
	echo 3 0.999977910 >>$@


clean : 
	rm -f gamma.svg gamma.data factorials.data error.svg error.data errorTable.data main.exe sfuns.dll

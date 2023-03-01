Gamma.svg : gamma.gpi gamma.data factorials.data Error.svg 
	gnuplot $<

gamma.data : main.exe
	mono $< -type:gamma > $@

Error.svg : error.gpi error.data errorTable.data
	gnuplot $<

error.data : main.exe
	mono $< -type:error > $@

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

errorTable.data : Makefile
	>$@
	echo 0 0 >>$@
	echo 1 0.842700793 >>$@
	echo 2 0.995322265 >>$@
	echo 3 0.999977910 >>$@


clean : 
	rm -f gamma.svg gamma.data factorials.data error.svg error.data errorTable.data main.exe sfuns.dll

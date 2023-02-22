TIME = time --portability --append --output $@

N = 1e9
Out.txt : main.exe pfor.exe
	>$@
	$(TIME) mono main.exe -terms:$N -threads:1 >> $@
	$(TIME) mono main.exe -terms:$N -threads:2 >> $@
	$(TIME) mono main.exe -terms:$N -threads:3 >> $@
	$(TIME) mono main.exe -terms:$N -threads:4 >> $@
	$(TIME) mono pfor.exe -terms:$N -threads:4 >> $@

main.exe : main.cs
	mcs $<

pfor.exe : pfor.cs
	mcs $<

.PHONEY:clean
clean:
	$(RM) -rf *.txt *.exe

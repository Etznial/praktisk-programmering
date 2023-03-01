SDKDIR=$(HOME)/.dotnet/sdk/6.0.102# or whatever your sdk is
CSCPATH=$(SDKDIR)/Roslyn/bincore/csc.dll
NETSTANDARD=$(SDKDIR)/ref/netstandard.dll
CONFIG=$(SDKDIR)/vstest.console.runtimeconfig.json
DOTNET=DOTNET_CLI_TELEMETRY_OPTOUT=1; dotnet
CSC=$(DOTNET) $(CSCPATH) -reference:$(NETSTANDARD)
RUN=$(DOTNET) exec --runtimeconfig $(CONFIG)

Out.txt:main.exe
	$(RUN) $< > $@ 2>log
	echo "\ntall matrix:" >> $@
	$(RUN) $< 4 3 >> $@ 2>>log

main.exe: main.cs matlib.dll
	$(CSC) $< -reference:$(lastword $^)
matlib.dll: ../matrix/matrix.cs ../matrix/vector.cs \
	givensqr.cs gsqr.cs
	$(CSC) $^ -target:library -out:$@

clean:
	$(RM) Out* *.dll *.exe log*

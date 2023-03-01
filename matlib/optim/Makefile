all: Out.txt
Out.txt:main.exe
	mono $< >$@ 2>Log

main.exe:main.cs swarm.dll matlib.dll
	mcs $< -r:$(call clist,$(filter %.dll,$^))

LIBDIR=..
LIBSRC=\
	$(LIBDIR)/matrix/vector.cs\
	$(LIBDIR)/matrix/matrix.cs
matlib.dll:$(LIBSRC)
	mcs $^ -target:library -out:$@

swarm.dll: swarm.cs matlib.dll
	mcs -target:library -out:$@ \
	$(addprefix -r:,$(filter %.dll,$^)) $(filter %.cs,$^)

clean:
	$(RM) *.exe *.dll [Oo]ut* [Ll]og*

comma:=,
empty:=
space:=$(empty) $(empty)
clist = $(subst $(space),$(comma),$(1))

LIBDIR=..
LIBSRC=\
	$(LIBDIR)/matrix/vector.cs\
	$(LIBDIR)/matrix/matrix.cs\
	$(LIBDIR)/QR/givensqr.cs

all:Out.txt
Out.txt:main.exe
	mono $< >$@ 2>Log

main.exe:main.cs roots.dll matlib.dll
	mcs -t:exe -o:$@ $(addprefix -r:,$(filter %.dll,$^)) $(filter %.cs,$^)

matlib.dll:$(LIBSRC)
	mcs $^ -target:library -out:$@

roots.dll:jacobian.cs newton.cs broyden.cs matlib.dll
	src="$(filter %.cs,$^)"; libs="$(addprefix -r:,$(filter %.dll,$^))";\
	mcs -target:library -out:$@ $$libs $$src

clean:
	$(RM) *.exe *.dll [Oo]out* [Ll]og*

comma:=,
empty:=
space:=$(empty) $(empty)
clist = $(subst $(space),$(comma),$(1))

Out.txt : main.exe
	mono main.exe > Out.txt

#sfuns.dll : sfuns.cs
#	mcs $^ -target:library -out:sfuns.dll

#main.exe : main.cs sfuns
#	mcs main.cs -reference:sfuns.dll -target:exe -out:main.exe

main.exe : main.cs
	mcs main.cs -target:exe -out:main.exe


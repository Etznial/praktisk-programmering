Out.txt : main.exe
	mono main.exe > Out.txt

#math.dll : main.cs 
	#mcs $^ -target:library -out:math.dll
	#mcs static_hello.cs static_world.cs -target:library -out:hello.dll

#main.exe : main.cs math.dll
#	mcs main.cs -reference:math.dll -target:exe -out:main.exe

sfuns.dll : sfuns.cs
	mcs $^ -target:library -out:sfuns.dll

main.exe : main.cs sfuns.dll
	mcs main.cs -reference:sfuns.dll -target:exe -out:main.exe

clean :
	rm -f Out.txt sfuns.dll main.exe

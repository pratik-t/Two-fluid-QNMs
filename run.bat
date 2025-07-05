make
@REM cls
Main.exe < input.txt

@REM f2py -c --opt="-march=native -flto -O3 -floop-nest-optimize" --fcompiler=gnu95 --compiler=mingw32 -m test test.f95
@REM cls

@REM f2py -m test -h test.pyd test.f95
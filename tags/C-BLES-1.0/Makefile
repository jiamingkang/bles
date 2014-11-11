products: BLES5

BLES5: BLES5.o EMatrix.o Numbering.o ABFG.o Solve.o Sens.o Levels.o Input.o Output.o
	gfortran -m32 -o BLES5 -lblas -llapack -lma57d -lmetdum -larpack -lla01\
		BLES5.o EMatrix.o Numbering.o ABFG.o \
		Solve.o Sens.o Levels.o Output.o Input.o

BLES5.o: BLES5.c
	gcc -m32 -c -o BLES5.o BLES5.c

EMatrix.o: EMatrix.c
	gcc -m32 -c -o EMatrix.o EMatrix.c

Numbering.o: Numbering.c
	gcc -m32 -c -o Numbering.o Numbering.c

ABFG.o: ABFG.c
	gcc -m32 -c -o ABFG.o ABFG.c

Solve.o : Solve.c
	gcc -m32 -c -o Solve.o Solve.c

Sens.o : Sens.c
	gcc -m32 -c -o Sens.o Sens.c

Levels.o : Levels.c
	gcc -m32 -c -o Levels.o Levels.c

Input.o : Input.c
	gcc -m32 -c -o Input.o Input.c

Output.o : Output.c
	gcc -m32 -c -o Output.o Output.c

clean: 
	rm *.o
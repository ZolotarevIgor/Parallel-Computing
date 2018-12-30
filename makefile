all:
	echo "use \"make ex***\" or name your file"

ex01a: ex01a.c
	gcc -o ex01a.px -O2 -lm ex01a.c
ex01b: ex01b.c
	gcc -o ex01b.px -O2 -lm ex01b.c
ex01c: ex01c.c
	gcc -o ex01c.px -O2 -lm ex01c.c
ex01d: ex01d.c
	gcc -o ex01d.px -O2 -lm ex01d.c
ex02a: ex02a.c mycom.c
	gcc -o ex02a.px -O2 -lm ex02a.c mycom.c
ex02b: ex02b.c mycom.c
	gcc -o ex02b.px -O2 -lm ex02b.c mycom.c
ex03а: ex03а.c
	gcc -o ex03a.px -O2 -lm ex03а.c
ex03b: ex03b.c
	gcc -o ex03b.px -O2 -lm ex03b.c
ex04a: mycom.c ex04a.c
	gcc -o ex04a.px -O2 -lm ex04a.c mycom.c
ex04b: mycom.c ex04b.c
	gcc -o ex04b.px -O2 -lm –pthread ex04b.c mycom.c
ex05а: mycom.c ex05а.c
	gcc -o ex05a.px -O2 -lm ex05а.c mycom.c
ex05b: mycom.c ex05b.c
	gcc -o ex05b.px -O2 -lm -pthread ex05b.c mycom.c
ex05c: mycom.c ex05c.c
	gcc -o ex05c.px -O2 -lm -fopenmp ex05c.c mycom.c

ex06a: ex06a.c
	mpicc -o ex06a.px -O2 -lm ex06a.c
	#mpirun -np 1 -nolocal -machinefile hosts ex06a.px
ex06b: ex06b.c
	mpicc -o ex06b.px -O2 -lm ex06b.c
 	#mpirun -np 1 -nolocal -machinefile hosts ex06b.px
ex06c: ex06c.c
	 mpicc -o ex06c.px -O2 -lm  ex06c.c
	 #mpirun -np 1 -nolocal -machinefile hosts ex06c.px
ex07a: ex07a.c mycom.c mynet.c
	 mpicc -o ex07a.px -O2 -lm  ex07a.c mycom.c mynet.c
	 #mpirun -np 1 -nolocal -machinefile hosts ex07a.px
ex08a: mynet.c ex08a.c
	mpicc -o ex08a.px -O2 -lm  ex08a.c mynet.c
	#mpirun -np 1 -nolocal -machinefile hosts ex08a.px 1 1 0
ex08b: mycom.c mynet.c ex08b.c
	mpicc -o ex08b.px -O2 -lm  ex08b.c mycom.c mynet.c
	#mpirun -np 1 -nolocal -machinefile hosts ex08b.px


ex09a: mynet.c ex09a.c
	mpicc -o ex09a.px -O2 -lm  ex09a.c mynet.c
 	#mpirun -np 2 -nolocal -machinefile hosts ex09a.px
ex09b: mynet.c ex09b.c
	 mpicc -o ex09b.px -O2 -lm  ex09b.c mynet.c
	 #mpirun -np 2 -nolocal -machinefile hosts ex09b.px
ex10a: mynet.c myrand.c ex10a.c
	 mpicc -o ex10a.px -O2 -lm  ex10a.c mynet.c myrand.c
	#mpirun -np 1 -nolocal -machinefile hosts ex10a.px
 ex11a: mycom.c mynet.c ex11a.c
	mpicc -o ex11a.px -O2 -lm  ex11a.c mycom.c mynet.c myprog.c
	#mpirun -np 1 -nolocal -machinefile hosts ex11a.px <10...10000000>
	#mpirun -np <1...16> -nolocal -machinefile hosts ex11a.px 10000000
 ex11b: mycom.c mynet.c  ex11b.c
	mpicc -o ex11b.px -O2 -lm  ex11b.c mycom.c mynet.c myprog.c
	#mpirun -np <1...16> -nolocal -machinefile hosts ex11a.px 10000000
 ex12a: mycom.c mynet.c  ex12a.c
	mpicc -o ex12a.px -O2 -lm  ex12a.c mycom.c mynet.c myprog.c
	#mpirun -np 1 -nolocal -machinefile hosts ex12a.px <10...1000000>
	#mpirun -np <1...6> -nolocal -machinefile hosts ex12a.px 1000000
 ex13a: mycom.c mynet.c myio.c  ex13a.c
	mpicc -o ex13a.px -O2 -lm  ex13a.c mycom.c mynet.c myio.c
	#mpirun -np 1 -nolocal -machinefile hosts ex13a.px
	#mpirun -np <1...12> -nolocal -machinefile hosts ex13a.px 10000 2000 1000
ex13b: mycom.c mynet.c ex13b.c
	 mpicc -o ex13b.px -O2 -lm  ex13b.c mycom.c mynet.c myio.c myprog.c
	 #mpirun -np 1 -nolocal -machinefile hosts ex13b.px
 ex14a: mycom.c mynet.c myio.c  ex14a.c
	mpicc -o ex14a.px -O2 -lm  ex14a.c mycom.c mynet.c myio.c
	#mpirun -np <1...16> -nolocal -machinefile hosts ex14a.px 100 100 5000 25000
ex14b: mycom.c mynet.c myio.c ex14b.c
	mpicc -o ex14b.px -O2 -lm  ex14b.c mycom.c mynet.c myio.c


clean:
	-rm *.p*; true
	-rm *.px; true

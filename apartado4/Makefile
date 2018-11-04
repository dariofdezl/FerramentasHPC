
all: dgesv

dgesv: dgesvfree.c
	gcc -g -p -lmkl -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -fopenmp -DREDUCTION -D ALIGNED -DNOALIAS -Ofast -march=native -o dgesv dgesvfree.c

run: 
	echo "Small test"
	./dgesv 2048
	echo "Medium test"
	./dgesv 4092
	echo "Large test"
	./dgesv 8192
clean: 
	rm dgesv

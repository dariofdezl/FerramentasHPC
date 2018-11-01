
all: dgesv

dgesv: dgesv.c
	icc -lmkl -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -g -p -fopenmp -Ofast -march=native -o dgesv dgesvfree.c

run: 
	echo "Small test"
	./dgesv 2048
	echo "Medium test"
	./dgesv 4096
	echo "Large test"
	./dgesv 8192
clean: 
	rm dgesv

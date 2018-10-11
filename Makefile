
all: dgesv

dgesv: dgesv.c
	icc -lmkl -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -g -p -fopenmp -O2  -march=native -o dgesv dgesv.c

run: 
	echo "Small test"
	./dgesv 1024

clean: 
	rm dgesv

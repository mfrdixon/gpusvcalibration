Error_Func.so: Error_Func.o SVFourierCosine.o

	nvcc -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -shared -o Error_Func.so Error_Func.o SVFourierCosine.o
Error_Func.o: Error_Func.c
	/opt/RevoR-6.2/bin/R CMD COMPILE Error_Func.c 
	gcc -std=gnu99 -I/opt/RevoR-6.2/lib64/Revo-6.2/R-2.15.3/lib64/R/include -DNDEBUG  -I/usr/local/include    -fpic  -g -O2  -c Error_Func.c -o Error_Func.o -fPIC
SVFourierCosine.o: SVFourierCosine.cu
	nvcc -gencode arch=compute_13,code=sm_13 -gencode arch=compute_20,code=sm_20 -c SVFourierCosine.cu -Xcompiler -fPIC

TESLA_SM13     := -gencode arch=compute_13,code=sm_13
FERMI_SM20     := -gencode arch=compute_20,code=sm_20
KEPLER_SM30    := -gencode arch=compute_30,code=sm_30 
KEPLER_SM35    := -gencode arch=compute_35,code=sm_35
ARCH_FLAGS     := $(TESLA_SM13) $(FERMI_SM20) $(KEPLER_SM30) $(KEPLER_SM35)
R_EXE          := /opt/RevoR-6.2/bin/R 
R_INCLUDE      := /opt/RevoR-6.2/lib64/Revo-6.2/R-2.15.3/lib64/R/include

Error_Func.so: Error_Func.o SVFourierCosine.o
	nvcc $(ARCH_FLAGS) -shared -o Error_Func.so Error_Func.o SVFourierCosine.o
Error_Func.o: Error_Func.c
	$(R_EXE) CMD COMPILE Error_Func.c
	gcc -std=gnu99 -I$(R_INCLUDE) -DNDEBUG  -I/usr/local/include    -fpic  -g -O2  -c Error_Func.c -o Error_Func.o -fPIC
SVFourierCosine.o: SVFourierCosine.cu		
	nvcc $(ARCH_FLAGS) -c SVFourierCosine.cu -Xcompiler -fPIC
clean :
	rm -f *~ *.o *.so

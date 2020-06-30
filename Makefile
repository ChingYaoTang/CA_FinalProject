TARGET	= MG_solver
CC	= g++ -fopenmp
CFLAGS	= -o
NVCC	= nvcc
OBJ	= main.o init.o basic.o prolongation.o restriction.o cal_residual.o relaxation.o relative_error.o exact_im.o up_down.o
RM	= rm -f

$(TARGET):$(OBJ)
	$(NVCC) $^ $(CFLAGS) $(TARGET)
main.o:main.cu
	$(NVCC) $(CFLAGS) $@ -c $^
init.o:init.c
	$(CC) $(CFLAGS) $@ -c $^
basic.o:basic.cu
	$(NVCC) $(CFLAGS) $@ -c $^
prolongation.o:prolongation.cu
	$(NVCC) $(CFLAGS) $@ -c $^
restriction.o:restriction.cu
	$(NVCC) $(CFLAGS) $@ -c $^
cal_residual.o:cal_residual.cu
	$(NVCC) $(CFLAGS) $@ -c $^
relaxation.o:relaxation.cu
	$(NVCC) $(CFLAGS) $@ -c $^
relative_error.o:relative_error.c
	$(CC) $(CFLAGS) $@ -c $^
exact_im.o:exact_im.c
	$(CC) $(CFLAGS) $@ -c $^
up_down.o:up_down.cu
	$(NVCC) $(CFLAGS) $@ -c $^



clean:
	$(RM) *.o $(TARGET)

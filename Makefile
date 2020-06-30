TARGET	= MG_solver
CC	= g++  -fopenmp
CFLAGS	= -o
OBJ	= main.o init.o basic.o prolongation.o restriction.o cal_residual.o relaxation.o relative_error.o exact_im.o up_down.o
RM	= rm -f

$(TARGET):$(OBJ)
	$(CC) $^ $(CFLAGS) $(TARGET)
main.o:main.c
	$(CC) $(CFLAGS) $@ -c $^
init.o:init.c
	$(CC) $(CFLAGS) $@ -c $^
basic.o:basic.c
	$(CC) $(CFLAGS) $@ -c $^
prolongation.o:prolongation.c
	$(CC) $(CFLAGS) $@ -c $^
restriction.o:restriction.c
	$(CC) $(CFLAGS) $@ -c $^
cal_residual.o:cal_residual.c
	$(CC) $(CFLAGS) $@ -c $^
relaxation.o:relaxation.c
	$(CC) $(CFLAGS) $@ -c $^
relative_error.o:relative_error.c
	$(CC) $(CFLAGS) $@ -c $^
exact_im.o:exact_im.c
	$(CC) $(CFLAGS) $@ -c $^
up_down.o:up_down.c
	$(CC) $(CFLAGS) $@ -c $^



clean:
	$(RM) *.o $(TARGET)

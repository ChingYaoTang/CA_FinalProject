TARGET	= MG_solver
CC	= mpic++
CFLAGS	= -o
OBJ	= main.o init_sin.o smoothing.o restriction.o basic.o prolongation.o cal_residual.o exactsolution.o relative_error.o
RM	= rm -f

$(TARGET):$(OBJ)
	$(CC) $^ $(CFLAGS) $(TARGET)
main.o:main.c
	$(CC) $(CFLAGS) $@ -c $^
init_sin.o:init_sin.c
	$(CC) $(CFLAGS) $@ -c $^
smoothing.o:smoothing.c
	$(CC) $(CFLAGS) $@ -c $^
restriction.o:restriction.c
	$(CC) $(CFLAGS) $@ -c $^
basic.o:basic.c
	$(CC) $(CFLAGS) $@ -c $^
prolongation.o:prolongation.c
	$(CC) $(CFLAGS) $@ -c $^
cal_residual.o:cal_residual.c
	$(CC) $(CFLAGS) $@ -c $^
exactsolution.o:exactsolution.c
	$(CC) $(CFLAGS) $@ -c $^
relative_error.o:relative_error.c
	$(CC) $(CFLAGS) $@ -c $^



clean:
	$(RM) *.o $(TARGET)

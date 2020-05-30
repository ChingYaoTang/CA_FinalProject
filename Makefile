TARGET	= MG_solver
CC	= mpic++
CFLAGS	= -o
OBJ	= main.o init_sin.o smoothing.o restriction.o basic.o prolongation.o
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


clean:
	$(RM) *.o $(TARGET)

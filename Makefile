TARGET	= MG_solver
CC	= mpic++
CFLAGS	= -o
OBJ	= main.o add.o init_sin.o
RM	= rm -f

$(TARGET):$(OBJ)
	$(CC) $^ $(CFLAGS) $(TARGET)
main.o:main.c
	$(CC) $(CFLAGS) $@ -c $^
add.o:add.c
	$(CC) $(CFLAGS) $@ -c $^
init_sin.o:init_sin.c
	$(CC) $(CFLAGS) $@ -c $^

clean:
	$(RM) *.o $(TARGET)

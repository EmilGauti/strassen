FLAGS = -Wall -O2
OBJS = main.o funcs.o
EXEC = matmul
CFLAGS=-Wall

all: $(EXEC)

$(EXEC): $(OBJS)
	gcc $(FLAGS) -o $(EXEC) $(OBJS)

funcs.o: funcs.h funcs.c
	gcc $(FLAGS) $(INCLUDES) -c funcs.c

main.o: main.c funcs.h
	gcc $(FLAGS) $(INCLUDES) -c main.c


clean:
	rm -f $(EXEC) $(OBJS)
CPP = /opt/homebrew/opt/llvm/bin/clang
CFLAGS = -O3 -Wall -Werror
CPPFLAGS = -I/opt/homebrew/opt/llvm/include -fopenmp
LDFLAGS = -L/opt/homebrew/opt/llvm/lib
RM = /bin/rm -f
OBJS = main.o funcs.o
EXECUTABLE = matmul

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS) $(CFLAGS)

funcs.o: funcs.h funcs.c
	$(CPP) $(CPPFLAGS) -c funcs.c $(CFLAGS)

main.o: main.c funcs.h
	$(CPP) $(CFLAGS) -c main.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS)

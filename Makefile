CC = gcc-9
CFLAGS = -Wall -Wextra -O3 #-O3 flag automatically slightly optimizes the code while compiling
LDFLAGS = -lm -lpng -pg -ggdb3  #Need an extra set of flags for functions.c (instead of just including -lm in CFLAGS)
SOURCES = functions.h functions.c main.c
OBJECTS = $(SOURCES:.c =.o)
EXECUTABLE = SysOLE

all : $(EXECUTABLE)

$(EXECUTABLE):$(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE) $(LDFLAGS) 

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
clean :
	rm -f $(EXECUTABLE)
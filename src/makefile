# définition des actions
.PHONY: clean

# définition des variables
CC = gcc
CCFLAGS = -O3 
LDFLAGS = 
SOURCES = spp_grasp.c
OBJECTS = $(SOURCES:.c =.o)

# regles
all: $(OBJECTS)
	$(CC) -o ../bin/spp_grasp $^ $(CCFLAGS) $(LDFLAGS)

%.o: %.c
	$(CC) -c $< $(CCFLAGS) $(LDFLAGS)

# action clean
clean:
	-rm *.o ../bin/spp_grasp *~

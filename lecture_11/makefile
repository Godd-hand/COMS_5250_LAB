CC = gcc
FFLAGS = -O3 
OBJECTS = main.o insert_term.o add.o evaluate_poly.o print.o delete_poly.o

.PHONY: clean

main.exe: $(OBJECTS)
	$(CC) $(OBJECTS) -o main.exe -lm

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe
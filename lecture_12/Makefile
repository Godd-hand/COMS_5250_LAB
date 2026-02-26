CC = gcc
FFLAGS = -O3 -Wall
OBJECTS = main.o bisector.o push.o pop.o peek.o DisplayStack.o GetStackSize.o DeleteStack.o

.PHONY: clean help

main.exe: $(OBJECTS)
	$(CC) $(FFLAGS) $(OBJECTS) -o main.exe -lm

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe

help:
	@echo "Valid targets:"
	@echo "  main.exe"
	@echo "  clean: removes *.o and .exe files"
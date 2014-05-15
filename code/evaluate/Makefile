.PHONY: all clean

all: evaluate

evaluate: evaluate.c
	gcc -Wall -lgmp -fopenmp -lgomp $< -o evaluate

clean:
	rm evaluate

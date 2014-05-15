.PHONY: all clean

all: evaluate

evaluate: evaluate.c
	g++ -Wall -lgmp -fopenmp -lgomp $< -o evaluate

clean:
	rm evaluate

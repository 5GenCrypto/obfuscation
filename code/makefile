CC 		= g++
CCFLAGS = -Wall --std=c++11 -g -O0 -fopenmp #-lrt
IFLAGS 	= -I/opt/local/include 
LFLAGS 	= -L/opt/local/lib -lgmpxx -lgmp

OBJS = src/clt_mlm.o src/utils.o 
SRCS = src/clt_mlm.cpp src/utils.cpp 

test.cpp: $(OBJS) $(SRCS) makefile
	$(CC) $(CCFLAGS) $(IFLAGS) $(LFLAGS) $(OBJS) test_mlm.cpp -o test_mlm

src/%.o: src/%.cpp
	$(CC) $(CCFLAGS) $(IFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS)
	rm -f test

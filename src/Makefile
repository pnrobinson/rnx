CC = g++
CCFLAGS = -g3 -ansi -DDEBUG 
LIBFLAGS = #lib/cppunitlite.a

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $<

objects = unittest.o Sequence.o Nussinov.o EnergyFunction2.o RNAStructure.o

all: maintest

.PHONEY= clean valgrind

cppunitlite:
	cd lib; make all

maintest: $(objects) maintest.o
	@echo "Compiling main unit test program..."
	${CC} -o maintest maintest.o $(objects) ${LIBFLAGS} 
	@echo "Running unit tests..."
	./maintest

valgrind:
	valgrind ./maintest

rnx: rnx.cpp $(objects)
	${CC} $(CCFLAGAS) -o rnx $(objects)    $<


clean:
	-rm $(objects) maintest maintest.o
	-rm *~
	-rm rnx

CC = g++
CCFLAGS = -g3 -ansi -DDEBUG 
LIBFLAGS = #lib/cppunitlite.a

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $<

objects = maintest.o unittest.o Sequence.o Nussinov.o EnergyFunction2.o RNAStructure.o

all: maintest

.PHONEY= clean valgrind

cppunitlite:
	cd lib; make all

maintest: $(objects)
	@echo "Compiling main unit test program..."
	${CC} -o maintest $(objects) ${LIBFLAGS} 
	@echo "Running unit tests..."
	#./maintest

valgrind:
	valgrind ./maintest

rnx: rnx.cpp
	${CC} $(CCFLAGAS) -o rnx $<


clean:
	-rm $(objects) maintest
	-rm *~
	-rm rnx

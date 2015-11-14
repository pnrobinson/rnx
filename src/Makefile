CC = g++
CCFLAGS = -g -ansi
LIBFLAGS = #lib/cppunitlite.a

%.o : %.cpp
	$(CC) $(CCFLAGS) -c $<

objects = maintest.o MyTestClass.o unittest.o Sequence.o

all: maintest

cppunitlite:
	cd lib; make all

maintest: $(objects)
	@echo "Compiling main unit test program..."
	${CC} -o maintest $(objects) ${LIBFLAGS} 
	@echo "Running unit tests..."
	./maintest
    
main.o: main.cpp
MyTestClass.o: MyTestClass.cpp
simplest_test.o: simplest_test.cpp
fixture_test.o: fixture_test.cpp


clean:
	rm $(objects) maintest
	-rm *~

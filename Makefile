

all: maintest

maintest:
	cd src; make maintest; ./maintest

clean:
	cd src; make clean
	-rm -rf rnx_docu

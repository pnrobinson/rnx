

all: maintest

maintest:
	cd src; make maintest; ./maintest

rnx:
	cd src; make rnx

clean:
	cd src; make clean
	-rm -rf rnx_docu

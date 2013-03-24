CC = gcc

#NUMPYINC = -I/Library/Frameworks/Python.framework/Versions/Current/include/python2.7
#ARRAYINC = -I/Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/numpy/core/include/numpy

NUMPYINC = -I/usr/include/python2.7
ARRAYINC = -I/opt/python/lib/python2.7/site-packages/numpy/core/include/numpy
#CFLAGS = -fpic -c $(NUMPYINC) $(ARRAYINC) -std=gnu99 -fnested-functions
CFLAGS = -c $(NUMPYINC) $(ARRAYINC) -std=gnu99 -fnested-functions
LIB = -L/usr/lib/python2.7/config -L/usr/local/lib

.SUFFIXES: .o .c
%.o: %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

default: MakePower.so

MakePower.so: cMakePower.o wrapper.o 
	$(CC) -shared -o $@ $^

clean:
	rm -f *.o *.c~

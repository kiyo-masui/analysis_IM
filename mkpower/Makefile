CC = gcc

NUMPYINC = -I/opt/python/2.7.3-mkl/include/python2.7/
ARRAYINC = -I/opt/python/2.7.3-mkl/lib/python2.7/site-packages/numpy/core/include/numpy
CFLAGS = -fpic -c $(NUMPYINC) $(ARRAYINC) -std=gnu99
LIB = -L/opt/python/2.7.3-mkl/lib/python2.7/config/ -L/usr/local/lib

.SUFFIXES: .o .c
%.o: %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

default: MakePower.so

MakePower.so: cMakePower.o wrapper.o 
	$(CC) -shared -o $@ $^

clean:
	rm -f *.o *.c~

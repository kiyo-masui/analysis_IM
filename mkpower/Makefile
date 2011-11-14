CC = gcc

#ARRAYINC = -I/usr/lib64/python2.6/site-packages/numpy/core/include/numpy/

NUMPYINC = -I/usr/include/python2.6/
ARRAYINC = -I/opt/numpy/lib64/python2.6/site-packages/numpy/core/include/numpy
CFLAGS = -fpic -c $(NUMPYINC) $(ARRAYINC) -std=gnu99
LIB = -L/usr/lib/python2.6/config -L/usr/local/lib

.SUFFIXES: .o .c
%.o: %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

default: MakePower.so

MakePower.so: cMakePower.o wrapper.o 
	$(CC) -shared -o $@ $^

clean:
	rm -f *.o *.c~

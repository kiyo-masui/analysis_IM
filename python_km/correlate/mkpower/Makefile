CC = gcc
# These need to be generalized using envornoment variables.
#NUMPYINC = -I/usr/include/python2.6/numpy
#CFLAGS = -fpic -c -I/usr/include/python2.6 $(NUMPYINC) -std=gnu99
NUMPYINC = -I/opt/sage/4.4.4/local/lib/python2.6/site-packages/numpy/core/include/numpy
CFLAGS = -fpic -c -I/opt/sage/4.4.4/local/include/python2.6 $(NUMPYINC) -std=gnu99
LIB = -L/usr/lib/python2.6/config -L/usr/local/lib

.SUFFIXES: .o .c
%.o: %.c
	$(CC) $(CFLAGS) $(LIB) -c $<

default: MakePower.so

MakePower.so: cMakePower.o wrapper.o 
	$(CC) -shared -o $@ $^

clean:
	rm -f *.o *.c~

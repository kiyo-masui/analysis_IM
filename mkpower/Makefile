CC = gcc

#ARRAYINC = -I/usr/lib64/python2.6/site-packages/numpy/core/include/numpy/

#NUMPYINC = -I/usr/include/python2.4/
#NUMPYINC = -I/obj/backup1/opt/python/2.7.1/include/python2.7
NUMPYINC = -I/opt/sage/4.4.4/local/include/python2.6
#ARRAYINC = -I/cita/h/home-2/xuelei/opt/numpy/lib/python2.6/site-packages/numpy/core/include/numpy/
ARRAYINC = -I/opt/sage/4.4.4/local/lib/python2.6/site-packages/numpy/core/include/numpy
#ARRAYINC = -I/opt/numpy/lib64/python2.6/site-packages/numpy/core/include/numpy
#ARRAYINC = -I/cita/h/home-3/kiyo/programs/python/numpy-1.5.0b1/numpy/core/include/numpy
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

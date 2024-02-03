CC = gcc
CFLAGS = -g -Wall -std=c99 -fopenmp -mavx -mfma -pthread
LDFLAGS = -fopenmp
# CUNIT = -L/home/ff/cs61c/cunit/install/lib -I/home/ff/cs61c/cunit/install/include -lcunit
# CUNIT = -L/usr/local/include/CUnit -I/usr/local/include/CUnit -lcunit
# CUNIT = -L/usr/local/lib -lcunit
CUNIT = -lcunit
# PYTHON = -I/usr/include/python3.6 -lpython3.6m
# PYTHON = -I/usr/local/include/python3.6dm -lpython3.6dm
PYTHON = -I/usr/include/python3.10 -lpython3.10

install:
	if [ ! -f files.txt ]; then touch files.txt; fi
	rm -rf build
	xargs rm -rf < files.txt
	python3 setup.py install --record files.txt

uninstall:
	if [ ! -f files.txt ]; then touch files.txt; fi
	rm -rf build
	xargs rm -rf < files.txt

clean:
	rm -f *.o
	rm -f test
	rm -rf build
	rm -rf __pycache__

test:
	rm -f test
	$(CC) $(CFLAGS) mat_test.c matrix.c -o test $(LDFLAGS) $(CUNIT) $(PYTHON)
	./test

# a:
# 	rm -f a 
# 	$(CC) $(CFLAGS) a.c -o a $(LDFLAGS) $(PYTHON)
# 	./a

.PHONY: test
CFLAGS += -Wall -W -O3 

default: mtest

mtest: mtest.o
	$(CC) mtest.o -ltommath -o mtest

clean:
	rm -f *.o mtest *~

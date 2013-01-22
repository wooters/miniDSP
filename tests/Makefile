
GP= $(HOME)/src/gnuplot_i/src
GP_OBJ=$(GP)/gnuplot_i.o
SD_SRCS = sinedriver.c 
SD_OBJS = $(SD_SRCS:.c=.o) $(GP_OBJ)

TLA_SRCS = testlocalaudio.c
TLA_OBJS = $(TLA_SRCS:.c=.o)

CFLAGS = -std=c99 -O3 -Wall

%.o:%.c
	$(CC) $(CFLAGS) -I $(GP) -I .. -c $*.c -o $@ 

default: sinedriver

sinedriver: $(SD_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SD_OBJS) $(LDFLAGS) -L.. -lminidsp -lfftw3

testlocalaudio: $(TLA_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TLA_OBJS) $(LDFLAGS) -L.. -lminidsp -lfftw3 -lportaudio

.PHONY:	clean
clean:
	-rm *.o sinedriver

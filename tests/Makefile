


GP= $(HOME)/src/gnuplot_i/src
GP_OBJ=$(GP)/gnuplot_i.o
SD_SRCS = sinedriver.c 
SD_OBJS = $(SD_SRCS:.c=.o) $(GP_OBJ)

CFLAGS = -std=c99 -O3 -Wall

%.o:%.c
	$(CC) $(CFLAGS) -I $(GP) -I .. -c $*.c -o $@ 

default: sinedriver

sinedriver: $(SD_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SD_OBJS) $(LDFLAGS) -L.. -lminidsp -lfftw3

.PHONY:	clean
clean:
	-rm *.o sinedriver

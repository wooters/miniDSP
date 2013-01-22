
GP= $(HOME)/src/gnuplot_i/src
GP_OBJ=$(GP)/gnuplot_i.o

GPT_SRCS = gcc_phat_test.c 
GPT_OBJS = $(GPT_SRCS:.c=.o) $(GP_OBJ)

GPMT_SRCS = gcc_phat_multi_test.c 
GPMT_OBJS = $(GPMT_SRCS:.c=.o)

TLA_SRCS = testliveio.c
TLA_OBJS = $(TLA_SRCS:.c=.o)

TLR_SRCS = testliverec.c $(GP_OBJ)
TLR_OBJS = $(TLR_SRCS:.c=.o)

CFLAGS = -O3 -Wall

%.o:%.c
	$(CC) $(CFLAGS) -I $(GP) -I .. -c $*.c -o $@ 

default: gcc_phat_test

gcc_phat_test: $(GPT_OBJS) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $(GPT_OBJS) $(LDFLAGS) -L.. -lminidsp -lfftw3

gcc_phat_multi_test: $(GPMT_OBJS) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $(GPMT_OBJS) $(LDFLAGS) -L.. -lminidsp -lfftw3

testliveio: $(TLA_OBJS) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $(TLA_OBJS) $(LDFLAGS) -L.. -lminidsp -lportaudio

testliverec: $(TLR_OBJS) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $(TLR_OBJS) $(LDFLAGS) -L.. -lminidsp -lportaudio

.PHONY:	clean
clean:
	-rm *.o gcc_phat_test gcc_phat_multi_test testliveio testliverec

CC = gcc
CFLAGS := -Wall -O3 -pipe -pg 

all: mapper mapper_sse mapper_pth
	date > compile_time

mapper: mapper.c iobank.h
	$(CC) -o mapper mapper.c $(CFLAGS)  -lm

mapper_sse: mapper_sse.c  iobank.h
	$(CC) -o mapper_sse mapper_sse.c $(CFLAGS)  -lm -msse4

mapper_pth: mapper_pth.c iobank.h
	$(CC) -o mapper_pth mapper_pth.c $(CFLAGS)  -lm -msse4 -lpthread -D_REENTRANT


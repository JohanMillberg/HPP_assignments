CC = gcc
LD = gcc
CFLAGS = -Wall
LDFLAGS =
LDLIBS =
RM = /bin/rm -f
OBJS = readData.o main.o sort.o
SORT = sorting
all:$(SORT)

$(SORT): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(SORT)

readData.o: readData.h readData.c
	$(CC) $(CFLAGS) -c readData.c

sort.o: sort.c sort.h
	$(CC) $(CFLAGS) -c sort.c

main.o: main.h main.c readData.h sort.h
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) $(SORT) $(OBJS)

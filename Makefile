CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		
PROG=		ropebwt bcrbwt
INCLUDES=	
LIBS=		-lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ropebwt:bprope6.o main.o
		$(CC) $(CFLAGS) $(DFLAGS) bprope6.o main.o -o $@ $(LIBS)

bcrbwt:bcr.o rld.o
		$(CC) $(CFLAGS) $(DFLAGS) bcr.o rld.o -o $@ $(LIBS)

rld.o:rld.c rld.h
		$(CC) -c $(CFLAGS) -D_NO_UTILS_H rld.c -o $@

bprope6.o:bprope6.h

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

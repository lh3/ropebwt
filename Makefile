CC=			gcc
CFLAGS=		-g -Wall -O2 #-fno-inline-functions -fno-inline-functions-called-once
DFLAGS=		
OBJS=		bprope6.o rld.o
PROG=		ropebwt
INCLUDES=	
LIBS=		-lpthread -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ropebwt:$(OBJS) ropebwt.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) ropebwt.o -o $@ $(LIBS)

rld.o:rld.c rld.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_DNA_ONLY -D_NO_UTILS_H rld.c -o $@

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

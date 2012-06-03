CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		
OBJS=		bprope6.o rbrope6.o rbrope6-mt.o
PROG=		ropebwt
INCLUDES=	
LIBS=		-lpthread -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

ropebwt:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

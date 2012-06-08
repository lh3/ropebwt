CC=			gcc
CFLAGS=		-g -Wall -O2 #-fno-inline-functions -fno-inline-functions-called-once
DFLAGS=
PROG=		bpr-mt ropebwt
INCLUDES=	
LIBS=		-lpthread -lz

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bpr-mt:bprope6.o rld.o bpr-mt.o main-bpr-mt.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

ropebwt:bprope6.o rbrope6.o bcr.o main-misc.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

bprope6.o:bprope6.h
rbrope6.o:rbrope6.h
bcr.o:bcr.h

rld.o:rld.c rld.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -D_DNA_ONLY -D_NO_UTILS_H rld.c -o $@

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

CC=		gcc

CFLAGS=		-Wall -Wextra -O3 -fomit-frame-pointer -ffast-math -std=gnu99 -march=native -mtune=native
DFLAGS=

#CFLAGS=		-Wall -Wextra -O0 -ggdb -std=gnu99
#DFLAGS=		-DPRGDEBUG

PROG=		semEP
SOLVER=		util.c types.c graph.c hash_map.c\
		CA.c metric.c semEP.c input.c main.c

SOLVEROBJS=	$(SOLVER:.c=.o)
INCLUDES=	-I.
LIBS=		-lm
INSTALLDIR=	../

.SUFFIXES:.c .o

all:		$(PROG)

$(PROG):	$(SOLVEROBJS)
		$(CC) $(CFLAGS) -o $(INSTALLDIR)$(PROG) $(SOLVEROBJS) $(LIBS)

.c.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) $(DFLAGS) $< -o $@

.PHONY : clean

clean :
	rm -rf $(INSTALLDIR)$(PROG) *.o *.dSYM *~

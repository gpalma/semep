CC=		cc

CFLAGS=		-Wall -Wextra -O3 -ffast-math -std=gnu99 -mtune=native
DFLAGS=

#CFLAGS=		-Wall -Wextra -O0 -ggdb -std=gnu99
#DFLAGS=		-DPRGDEBUG

PROG=		semEP

SOLVER=		util.c rbtree.c graph_adj.c hash_map.c input.c semEP.c main.c
LIBS=		-lm
SOLVEROBJS=	$(SOLVER:.c=.o)

INSTALLDIR=	../

.SUFFIXES:.c .o

all:		$(PROG)

$(PROG):	$(SOLVEROBJS)
		$(CC) $(CFLAGS) -o $(INSTALLDIR)$(PROG) $(SOLVEROBJS) $(LIBS)

.c.o:
	$(CC) -c  $(CFLAGS) $(DFLAGS) $< -o $@

.PHONY : clean

clean :
	rm -rf $(INSTALLDIR)$(PROG) *.o *.dSYM $(INSTALLDIR)*.dSYM *~

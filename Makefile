# read npc from ../generatempi.f line n.1
NPC := $(shell head -n 1 ../generatempi.f | awk -F 'npc=' '{print $$2}' | awk -F ',' '{print $$1}')
CFLAGS := -O3 -D NPC=$(NPC)

all: lista isf struct msd fself

# and build the list of file,index,time for each configuration ../Cnf*-idx
lista:
	bash lista.sh

# compila isf.c -> isf.x
isf:
	gcc ${CFLAGS} -D MESHFACT=1 isf.c -lm -o isf.x

struct: 
	gcc ${CFLAGS} -D MESHFACT=1 struct.c -lm -o struct.x

msd:
	gcc ${CFLAGS} msd-lin.c -lm -o msd-lin.x
	gcc ${CFLAGS} msd-log.c -lm -o msd-log.x

fself:
	gcc ${CFLAGS} -D MESHFACT=1 fself.c -lm -o fself.x

wavenumbers:
	gcc -O3 -D MESHFACT=1 wavenumbers.c -lm -o wavenumbers.x -Wall

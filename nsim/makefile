OBJECTS  = nsim.o sim.o net.o mob.o rl.o stat.o data.o heap.o
.SUFFIXES : .o .cpp

CC = g++

OPTION = -g

CFLAGS  = -c $(OPTION)

nsim : $(OBJECTS)
	$(CC)  -o nsim $(OBJECTS)

sim.o : sim.cpp sim.h global.h net.h mob.h rl.h data.h
nsim.o : nsim.cpp sim.h
net.o : net.cpp net.h mob.h rl.h
mob.o : mob.cpp mob.h net.h rl.h stat.h
rl.o : rl.cpp rl.h net.h ef.h
stat.o : stat.cpp stat.h list.h net.h mob.h
data.o : data.cpp data.h net.h list.h mob.h
heap.o : heap.cpp heap.h list.h

.cpp.o :
	$(CC) $(CFLAGS) $<

#include "heap.h"
#include <memory.h>
#include "sim.h"




Heap::Heap(int max) : c(max), last(-1)
{
	create(&ev, c);
}


void Heap::put(Event* event)
{
	int i = 1 + last++;
	if (i == c) {
		ev = grow(ev, c, INC);
		c +=INC;
	}
	int par = (i - 1)/2;	// parent
	while (i > 0 && event->t < ev[par]->t) {
		ev[i] = ev[par];
		i = par;
		par = (i - 1)/2;
	}
	ev[i] = event;
}


Event* Heap::get()
{
	Event* first = *ev;
	Event* lastEv = ev[last--];
	/* adjust heap by swapping earlier of children with 
	 parent until correct place for bottom on heap is found */
	int par = 0;
	int child = earlier(par);
	while (child && lastEv->t > ev[child]->t) {
		ev[par] = ev[child];
		par = child;
		child = earlier(par);
	}
        if (lastEv != first) 
	  ev[par] = lastEv;
	else
	  ev[par] = 0;
	  
	return first;
}


int Heap::earlier(int par) {
	int lchild = 2*par + 1;
	if (lchild > last)
		return 0;
	int rchild = lchild + 1;
	if (rchild <= last && ev[lchild]->t > ev[rchild]->t)
		return rchild;
	else
		return lchild;
}


/* void Heap::grow(int add)
{
	int newC = c + add;
	Event* newEv = Event*[newC];
	memset(newEv,0,newC*sizeof(Event);
	memcpy(newEv,ev,c*sizeof(Event*));
	delete ev;
	ev = newEv;
	c = newC;
}	 */

	

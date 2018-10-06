#ifndef HEAP_H
#define HEAP_H


const int INC = 100;


struct Event {
	double t;
	int	type;
	int	id;
	struct Mob*	gen;  // event generator
	Event(double _t, int _type, Mob* _g, int _id = 0) :
		t(_t), type(_type), gen(_g), id(_id) {};
};

struct Heap {
	Event** ev;
	int last;
	int c;
	void put(Event*);
	Event* get();
	int earlier(int);
	Heap(int max);
	
};



template <class Type>
	Type* grow(Type *old, int s, int add)
{
	int s1 = s + add;
	Type* newAr = new Type[s1];
	memset(newAr, 0, s1*sizeof(Type));
	memcpy(newAr, old, s*sizeof(Type));
	delete old;
	return newAr;
}
	
// used for events and mobs in Net


#endif

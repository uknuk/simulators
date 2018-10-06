#include "list.h"
#include <stdio.h>


template <class Type> 
	void List<Type>::put(Node<Type>* node)
{
	if (!head)
		 head = node;
	else {
		tail->next = node;
		node->prev = tail;
	}
	tail = node;
	n++;
}


template <class Type>
	void List::remove(Node<Type>* node)
{
	
	if (node == head) 
		head = node->next;
   if (node == tail) 
		tail = node->prev;
	if (node->prev)
			node->prev->next = node->next;
	if (node->next)
		node->next->prev = node->prev;
	n--;
	
}



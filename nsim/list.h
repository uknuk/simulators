#ifndef LIST_H
#define LIST_H



template <class Type> struct Node {
	Type item;
	Node* next;
	Node* prev;
	Node(Type i) : item(i), next(0), prev(0) {};
}; 

template <class Type> struct List{
	int n;
	Node<Type>* head;
	Node<Type>* tail;
	List() : n(0) , head(0), tail(0) {};
	void put(Node<Type>*);
	void remove(Node<Type>*);
};

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
	void List<Type>::remove(Node<Type>* node)
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


#endif




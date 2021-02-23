// Method Definitions for ILLNode

template <class T>
void ILLNode<T>::remove() {
    prev->next = next;
    next->prev = prev;
    next = this;
    prev = this;
};

template <class T>
void ILLNode<T>::insert_after(ILLNode<T> &node) {
    remove();
    prev = &node;
    next = node.next;
    node.next = this;
    next->prev = this;
};
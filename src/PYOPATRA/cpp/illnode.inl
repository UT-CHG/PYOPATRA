// Method Definitions for ILLNode

#include <iostream>

template <class T>
void ILLNode<T>::remove() {
    if (prev != nullptr) {
        prev->next = next;
    }

    if (next != nullptr) {
        next->prev = prev;
    }

    next = nullptr;
    prev = nullptr;
}

template <class T>
void ILLNode<T>::insert_after(ILLNode<T> &node) {
    remove();
    prev = &node;
    next = node.next;
    node.next = this;

    if (next != nullptr) {
        next->prev = this;
    }
}
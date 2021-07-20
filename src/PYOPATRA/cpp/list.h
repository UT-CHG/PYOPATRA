//
// Created by Georgia Stuart on 2/23/21.
//

#ifndef PYOPATRA_LIST_H
#define PYOPATRA_LIST_H

#include "illnode.h"

template <class T>
class List {
public:
    size_t length{};

    List() : length(0), head(nullptr), tail(nullptr){};


    void push(ILLNode<T> &node);
    void remove(ILLNode<T> &node);
//    ILLNode<T>* pop_current();
//    void advance();
    ILLNode<T>* get_head() { return head; }
    ILLNode<T>* get_tail() { return tail; }
    const ILLNode<T>* get_head() const { return head; }
    const ILLNode<T>* get_tail() const { return tail; }

private:
    ILLNode<T> *head, *tail;

};

#include "list.inl"

#endif //PYOPATRA_LIST_H

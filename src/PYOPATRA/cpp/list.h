//
// Created by Georgia Stuart on 2/23/21.
//

#ifndef PYOPATRA_LIST_H
#define PYOPATRA_LIST_H

#include "illnode.h"

template <class T>
class List {
public:
    ILLNode<T> *head, *tail, *current;
    size_t length{};

    List() : head(nullptr), tail(nullptr), current(nullptr), length(0) {};


    void push(ILLNode<T> &node);
    void remove(ILLNode<T> &node);
    ILLNode<T>* pop_current();
    void advance();

};

#include "list.inl"

#endif //PYOPATRA_LIST_H

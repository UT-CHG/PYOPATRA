//
// Created by Georgia Stuart on 2/23/21.
//

#ifndef PYOPATRA_LIST_H
#define PYOPATRA_LIST_H

#include "illnode.h"

template <class T>
class List {
public:
    ILLNode<T> *head, *tail
    size_t length;

    List() : head(nullptr), tail(nullptr), length(0) {};
};

#endif //PYOPATRA_LIST_H

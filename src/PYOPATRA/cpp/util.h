//
// Created by georgia on 7/26/21.
//

#ifndef PYOPATRA_UTIL_H
#define PYOPATRA_UTIL_H

#include <random>

static std::normal_distribution<double> normal(0, 1);
static std::uniform_real_distribution<double> unif_pi(0, M_PI);
static std::uniform_real_distribution<double> unif(0, 1);
static std::default_random_engine generator;

template <typename T>
class PointerWrapper {
private:
    T *ptr;

public:
    PointerWrapper() = default;
    PointerWrapper(T* ptr)
        : ptr(ptr)
    {}

    T* get_pointer() { return ptr; }
    void set_pointer(T* new_ptr) { ptr = new_ptr; }
};

#endif //PYOPATRA_UTIL_H

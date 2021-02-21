//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef LPTCPP_COORDINATE_H
#define LPTCPP_COORDINATE_H

template<typename T>
class CoordinateT {
public:
    T x, y, z;

    CoordinateT() : x(0), y(0), z(0) { }
    CoordinateT(T x, T y, T z) : x(x), y(y), z(z) { }
};


typedef CoordinateT<double> CoordinateD;
typedef CoordinateT<float> CoordinateF;
typedef CoordinateT<int> CoordinateI;
#endif //LPTCPP_COORDINATE_H

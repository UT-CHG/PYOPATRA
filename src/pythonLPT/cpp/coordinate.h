//
// Created by Georgia Stuart on 2/4/21.
//

#ifndef LPTCPP_COORDINATE_H
#define LPTCPP_COORDINATE_H

template<typename T>
class CoordinateT {
public:
    T latitude, longitude, depth;

    CoordinateT() : latitude(0), longitude(0), depth(0) { }
    CoordinateT(T latitude, T longitude, T depth) : latitude(latitude), longitude(longitude), depth(depth) { }
};


typedef CoordinateT<double> CoordinateD;
typedef CoordinateT<float> CoordinateF;
typedef CoordinateT<int> CoordinateI;
#endif //LPTCPP_COORDINATE_H

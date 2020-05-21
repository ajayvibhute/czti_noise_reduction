/* 
 * @file Mvector.h
 * @author Tanul Gupta
 * @date Created on Oct 14, 2015, 9:20 AM
 * @brief Class to operate with maths vectors
 * @details 
 */

#ifndef MVECTOR_H
#define	MVECTOR_H

#include <iostream>
#include <cmath>
#include <vector>
#include "macrodef.h"
#include "glog/logging.h"

using namespace std;

//Maths vector class
class MVector {
private:
    double x;
    double y;
    double z;
    double thetaRad;
    double thetaDeg;
    double phiRad; //measured from the z axis
    double phiDeg;
    double phiDashRad; //measured from the xy plane.
    double phiDashDeg;
    double magnitude;

public:
    MVector(double x, double y, double z);
    MVector(double thetaDeg, double phiDeg, bool measuredFromZ = false, double magnitude=1.0);

    double get_magnitude() {return magnitude;}
    double get_theta() {return thetaDeg;}
    double get_phi() {return phiDeg;}
    friend double dot_product(MVector a, MVector b);

    friend double get_angle_bw_vectors(MVector a, MVector b);
    void display();
};

#endif /**MVECTOR_H**/
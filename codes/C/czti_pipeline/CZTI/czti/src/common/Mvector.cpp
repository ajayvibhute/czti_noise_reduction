#include "Mvector.h"


//Maths vector class

MVector::MVector(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;

    //calculating magnitude
    this->magnitude = sqrt(x*x + y*y + z*z);
    //calculating theta and phi and storing them in class
    this->thetaRad = atan(y/x);
    this->thetaDeg = this->thetaRad * TODEG;
    
    this->phiRad = acos(z / this->magnitude);
    this->phiDeg = this->phiRad;
    this->phiDashDeg = 90 - phiDeg;
    this->phiDashRad = phiDashDeg*TORAD;
}

MVector::MVector(double thetaDeg, double phiDeg, bool measuredFromZ, double magnitude) {
    this->thetaDeg = thetaDeg;
    this->thetaRad = thetaDeg*TORAD;
    if(measuredFromZ==true){
        this->phiDeg = phiDeg;
        this->phiDashDeg = 90 - phiDeg;
        this->phiRad = phiDeg*TORAD;
        this->phiDashRad = phiDashDeg*TORAD;
    } else if (measuredFromZ==false){
        this->phiDashDeg = phiDeg;
        this->phiDeg = 90 - phiDeg;
        this->phiRad = phiDeg*TORAD;
        this->phiDashRad = phiDashDeg*TORAD;
    }
    this->magnitude = 1.0;
    //calculating x, y & z
    this->x = magnitude*(cos(thetaRad)*cos(phiDashRad));
    this->y = magnitude*(sin(thetaRad)*cos(phiDashRad));
    this->z = magnitude*sin(phiDashRad);
}

void MVector::display(){
    LOG(INFO) << "x: " << x;
    LOG(INFO) << "y: " << y;
    LOG(INFO) << "z: " << z;
    LOG(INFO) << "thetaDeg: " << thetaDeg;
    LOG(INFO) << "phiDeg: " << phiDashDeg;
}
double dot_product(MVector a, MVector b){
    double dotProductValue=0.0;
    dotProductValue = a.x * b.x + a.y*b.y + a.z*b.z;
    return dotProductValue;
}

double get_angle_bw_vectors(MVector a, MVector b){
    double angleRad = acos(dot_product(a,b)/(a.magnitude*b.magnitude));
    double angleDeg = angleRad*TODEG;
    return angleDeg;
}



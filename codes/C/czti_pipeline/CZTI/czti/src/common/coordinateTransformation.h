
/* 
 * @file coordinateTransformation.h
 * @author Tanul Gupta
 * @date May 23, 2015, 2:48 PM
 * @brief This file provides methods to transform coordinates from one system to another.
 * @details 
 */

#ifndef COORDINATETRANSFORMATION_H
#define	COORDINATETRANSFORMATION_H

#include "level1handler.h"
#include "caldbHandler.h"
#include "level2validation.h"
#include <vector>
#include <numeric>
#include "glog/logging.h"
#include "mkfRegeneration.h"

/**
 * Function to transform CZTI body vector into inertial vector
 * @param mkf
 * @param teldef
 * @param time
 * @param detX
 * @param detY
 * @param detZ
 * @param nX
 * @param nY
 * @param nZ
 * @return 
 */
int get_inertial_vector(Mkf &mkf, Teldef teldef, double time, 
                        float detX, float detY, float detZ, float &nX, float &nY, float &nZ);

/**
 * Function to calculate RA, DEC & TWIST values by taking average of 
 * inertial pointing and twist vectors.
 * @param vec_nX
 * @param vec_nY
 * @param vec_nZ
 * @param vec_nXt
 * @param vec_nYt
 * @param vec_nZt
 * @param RA
 * @param DEC
 * @param TWIST
 * @param status
 * @return 
 */
int calculate_RA_DEC_Twist(vector <float> vec_nX, vector <float> vec_nY, vector <float> vec_nZ
                            , vector <float> vec_nXt, vector <float> vec_nYt, vector <float> vec_nZt, 
                            float &RA, float &DEC, float &TWIST, int &status);

/**
 * @brief: Calculates RA, DEC & TWIST value for a single pointing and twist vector.
 * @param nX
 * @param nY
 * @param nZ
 * @param nXt
 * @param nYt
 * @param nZt
 * @param RA : in radians
 * @param DEC : in radians
 * @param TWIST : in radians
 * @param status
 * @return 
 */
int calculate_RA_DEC_Twist(float nX, float nY, float nZ, float nXt, float nYt, float nZt,
                            float &RA, float &DEC, float &TWIST, int &status);
/**
 * Function to transform inertial vector into CZTI body vector
 * @param mkf
 * @param teldef
 * @param nX
 * @param nY
 * @param nZ
 * @param detX
 * @param detY
 * @param detZ
 * @return 
 */
int get_body_vector(Mkf mkf, Teldef teldef, float nX, float nY, float nZ, float &detX, float &detY,
                    float &detZ);


/**
 * Function to convert detX, detY & detZ into thetaX and thetaY.
 * @param detX
 * @param detY
 * @param detZ
 * @param thetaX
 * @param thetaY
 * @return 
 */
int to_thetaX_thetaY(float detX, float detY, float detZ, float &thetaX, float &thetaY);


/**
 * Function to convert thetaX & thetaY into body vector (detX, detY & detZ)
 * @param thetaX
 * @param thetaY
 * @param detX
 * @param detY
 * @param detZ
 * @return 
 */
int to_detX_detY_detZ(float thetaX, float thetaY, float &detX, float &detY, float &detZ);

/**
 * Function to get nX, nY, nZ (inertial vector) values from RA and DEC (in radians).
 * @param RA
 * @param DEC
 * @param nX
 * @param nY
 * @param nZ
 * @return 
 */
int to_nX_nY_nZ(float RA, float DEC, float &nX, float&nY, float &nZ);

/**
 * Function to convert nX, nY, nZ (inertial vector) to RA & DEC (in radians).
 * @param nX
 * @param nY
 * @param nZ
 * @param RA
 * @param DEC
 * @return 
 */
int to_RA_DEC(float nX, float nY, float nZ, float &RA, float &DEC);
/**
 * Function to convert RA and DEC of source into Camera coordinates ThetaX and ThetaY. This function
 * only requires RA, DEC, TWIST (of pointing vector) to perform this operation.
 * @param aspectFilename
 * @param RAsource: radians
 * @param DECsource: radians
 * @param thetaX: radians
 * @param thetaY: radians
 * @return 
 */
template <class T> int to_thetaX_thetaY(string aspectFilename, T RAsource, T DECsource,
        T &thetaX, T &thetaY) {

    AspectFileHandler aspect;
    int status = 0; //status variable
    aspect.read_aspect_file(aspectFilename);
    double vX, vY, vZ = 0.0;
    double RA, DEC, TWIST = 0.0; // RA, DEC & TWIST of pointing vector to be stored in radians
    RA = aspect.get_RA() * M_PI / 180;
    DEC = aspect.get_DEC() * M_PI / 180;
    TWIST = aspect.get_TWIST() * M_PI / 180;
    DLOG(INFO) << "Read from aspect file RA:" << RA << "  DEC:" << DEC << " TWIST:" << TWIST;
    vZ = cos(DEC) * cos(DECsource) * cos(RAsource - RA) + sin(DEC) * sin(DECsource);
    vX = -cos(DECsource) * sin(RAsource - RA);
    vY = -(sin(DEC) * cos(DECsource) * cos(RAsource - RA) - sin(DECsource) * cos(DEC));

    thetaX = atan((vX * cos(TWIST) + vY * sin(TWIST)) / vZ);
    thetaY = atan((vY * cos(TWIST) - vX * sin(TWIST)) / vZ);
    return status;

}
/**
 * Function to convert RA & DEC of source into camera coordinates thetaX and thetaY.
 * @param RA in radians
 * @param DEC: in radians
 * @param TWIST: in radians
 * @param RAsource: in radians
 * @param DECsource: in radians
 * @param thetaX: in radians
 * @param thetaY: in radians
 * @return 
 */
template <class T>
int to_thetaX_thetaY(T RA, T DEC, T TWIST, T RAsource, T DECsource, 
        T &thetaX, T &thetaY){
    int status=0; //status variable
    double vX, vY, vZ = 0.0;
    

    vZ = cos(DEC)*cos(DECsource)*cos(RAsource-RA) + sin(DEC)*sin(DECsource);
    vX = -cos(DECsource)*sin(RAsource-RA);
    vY = -( sin(DEC)*cos(DECsource)*cos(RAsource-RA) - sin(DECsource)*cos(DEC));
    
    thetaX = atan((vX*cos(TWIST) + vY*sin(TWIST))/vZ);
    thetaY = atan((vY*cos(TWIST) - vX*sin(TWIST))/vZ);
    
    
    return status;
}

/**
 * Function to convert ThetaX, ThetaY into Celestial coordinates RA & DEC using pointing vector direction
 * from Aspect file.
 * @param aspectFileName
 * @param thetaX: in radians
 * @param thetaY: in radians
 * @param RAobject: in radians
 * @param DECobject: in radians
 * @return 
 */
template < class T> int to_RA_DEC(string aspectFileName, T thetaX, T thetaY, T &RAobject, T &DECobject){
    AspectFileHandler aspect;
    int status=0;
    double vZ, vX, vY=0.0;
    double xM, yM=0.0;
    double zS, xS, yS=0.0;
    double RA, DEC, TWIST=0.0;
    aspect.read_aspect_file(aspectFileName);

    RA = aspect.get_RA()* M_PI / 180;
    DEC = aspect.get_DEC()* M_PI / 180;
    TWIST = aspect.get_TWIST()* M_PI / 180;
    DLOG(INFO) << "Read from aspect file RA:" << RA << "  DEC:" << DEC << " TWIST:" << TWIST;
    vZ = 1/sqrt(1+ tan(thetaX)*tan(thetaX) + tan(thetaY)*tan(thetaY));
    vX= vZ*tan(thetaX);
    vY= vZ*tan(thetaY);
    
    xM=vX*cos(TWIST) -vY*sin(TWIST);
    yM=vY*cos(TWIST) + vX*sin(TWIST);
    
    zS = vZ*sin(DEC) + yM*cos(DEC);
    xS= vZ*cos(RA)*cos(DEC) - yM*cos(RA)*sin(DEC) + xM*sin(RA);
    yS = vZ*sin(RA)*cos(DEC) - yM*sin(RA)*sin(DEC) - xM*cos(RA);
    
    DECobject = asin(zS);

    if( atan2(yS, xS) < 0) {
        RAobject = atan2(yS, xS) + 2*M_PI;
    }
    else {
        RAobject = atan2(yS, xS);
    }    


    return status;
    
    
}

template <class T> 
int get_rpy_RA_DEC(Q q, T *rollRA, T* rollDEC, T* pitchRA, T* pitchDEC,
                                        T* yawRA, T* yawDEC){
    int status=0;
    int i=0; 
    Q quat;
    Axis axis, rot_axis;
    T RQ[3][3];
    T transRQ[3][3];
    T axis31[3][1];
    T sum=0.0;
    T rotatedAxis[3][1];
    T tempRA=0.0;
    T q1=0.0, q2=0.0, q3=0.0, q4=0.0;
    quat = q;

    q1 = quat.q1;
    q2 = quat.q2;
    q3 = quat.q3;
    q4 = quat.q4; //scalar
    
    //
    RQ[0][0] = q1 * q1 - q2 * q2 - q3 * q3 + q4*q4;
    RQ[0][1] = 2 * (q1 * q2 + q3 * q4);
    RQ[0][2] = 2 * (q1 * q3 - q2 * q4);
    RQ[1][0] = 2 * (q1 * q2 - q3 * q4);
    RQ[1][1] = -q1 * q1 + q2 * q2 - q3 * q3 + q4*q4;
    RQ[1][2] = 2 * (q2 * q3 + q1 * q4);
    RQ[2][0] = 2 * (q1 * q3 + q2 * q4);
    RQ[2][1] = 2 * (q2 * q3 - q1 * q4);
    RQ[2][2] = -q1 * q1 - q2 * q2 + q3 * q3 + q4*q4;

    //taking transpose
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            transRQ[i][j] = RQ[j][i];
        }
    }
    //transpose taken
    
    //calculating rollRA, rollDEC
    axis.x=0;
    axis.y=1;
    axis.z=0;
    axis31[0][0] = 0;
    axis31[1][0] = 1;
    axis31[2][0] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
            sum = 0;
            for (int k = 0; k < 3; k++) {
                sum = sum + (transRQ[i][k] * axis31[k][j]);
                rotatedAxis[i][j] = sum;
            }
        }
    }

    *rollDEC = (asin(rotatedAxis[2][0]))*180 / M_PI;
    tempRA = atan2(rotatedAxis[1][0], rotatedAxis[0][0]);
    tempRA = tempRA * 180 / M_PI;
    if (tempRA < 0) {
        *rollRA = tempRA + 360;
    } else {
        *rollRA = tempRA;
    }
    
    //calculating pitchRA, pitchDEC
    axis.x = 0;
    axis.y = 0;
    axis.z = 1;
    axis31[0][0] = 0;
    axis31[1][0] = 0;
    axis31[2][0] = 1;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
            sum = 0;
            for (int k = 0; k < 3; k++) {
                sum = sum + (transRQ[i][k] * axis31[k][j]);
                rotatedAxis[i][j] = sum;
            }
        }
    }

    *pitchDEC = (asin(rotatedAxis[2][0]))*180 / M_PI;
    tempRA = atan2(rotatedAxis[1][0], rotatedAxis[0][0]);
    tempRA = tempRA * 180 / M_PI;
    if (tempRA < 0) {
        *pitchRA = tempRA + 360;
    } else {
        *pitchRA = tempRA;
    }
    
    //calculating yawRA, yawDEC
    axis.x = 1;
    axis.y = 0;
    axis.z = 0;
    axis31[0][0] = 1;
    axis31[1][0] = 0;
    axis31[2][0] = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
            sum = 0;
            for (int k = 0; k < 3; k++) {
                sum = sum + (transRQ[i][k] * axis31[k][j]);
                rotatedAxis[i][j] = sum;
            }
        }
    }

    *yawDEC = (asin(rotatedAxis[2][0]))*180 / M_PI;
    tempRA = atan2(rotatedAxis[1][0], rotatedAxis[0][0]);
    tempRA = tempRA * 180 / M_PI;
    if (tempRA < 0) {
        *yawRA = tempRA + 360;
    } else {
        *yawRA = tempRA;
    }
    
    return EXIT_SUCCESS;
}
#endif	/* COORDINATETRANSFORMATION_H */


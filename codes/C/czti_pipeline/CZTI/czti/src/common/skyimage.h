/* 
 * @file  skyimage.h
 * @author Tanul Gupta
 * @date Created on May 8, 2015, 1:03 PM
 * @brief To create sky image from DPI and detect peaks in it
 * @details 
 * @version 0.2.0
 */

#ifndef SKYIMAGE_H
#define SKYIMAGE_H

#include "utils.h"

class Skyimage{
public:
    //Variables
    vector <long> dpi;
    vector <unsigned char> mask;
    void create_sky_image(vector <long> dpi, vector <unsigned char> mask);
};

//INDEPENDENT FUNCTIONS

#endif /*SKYIMAGE_H*/
/* 
 * File:   jpeg_handling.h
 * Author: Tanul Gupta
 *
 * Modified on November 19, 2014
 */
#ifndef JPEGHANDLING_H
#define JPEGHANDLING_H

#include <iostream>
#include "glog/logging.h"
#include <jpeglib.h>
#include <cstdlib>
#include <vector>
#include <fstream>

 using namespace std;

//>>>>>>>>>>>>>>>>>>> JPEG HANDLING

int makeJpeg(vector <unsigned char> &image, int height, int width, int input_components, J_COLOR_SPACE color_space, int quality, char* output_filename);
int makeJpeg(unsigned char* image, int height, int width, int input_components, J_COLOR_SPACE color_space, int quality, char* output_filename);

template<class T> int scaleDownToUC(T *TImage, unsigned char *UCimage, int height, int width){
    LOG(INFO)<< "-- Scaling down image to unsigned char --";
    int i=0, j=0; //counter variables;
    T temp_max=0;
    for(i=0; i<height*width; i++){
        if(TImage[i]>temp_max) { 
            temp_max=TImage[i]; 
        }
    }
    LOG(INFO)<<"MAXIMUM VALUE IS FOUND TO BE :"<< temp_max;
    for(i=0; i<height*width; i++){
        UCimage[i] = TImage[i]* 255 / temp_max;
    }
    
    LOG(INFO) << "--Scaling down completed --";
    
    return EXIT_SUCCESS;
}

//END >>>>>>>>>>>>>>JPEG HANDLING

#endif /*JPEGHANDLING_H*/
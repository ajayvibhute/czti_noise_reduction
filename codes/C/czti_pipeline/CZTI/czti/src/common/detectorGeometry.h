/* 
 * @file  detectorGeometry.h
 * @author Tanul Gupta
 * @date Created on June 25, 2015, 3:44 PM
 * @brief This class calculates coordinates, area, pixel locations etc. corresponding to CZTI detector
 * @details This class calculates coordinates, area, pixel locations etc. corresponding to CZTI detector
 * @version 1.0
 */

#ifndef DETECTORGEOMETRY_H
#define	DETECTORGEOMETRY_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include "utils.h"
#include "glog/logging.h"
#include "fitsio.h"
#include "macrodef.h"


class CztDetectorGeometry {
private:
    vector <float> xdet0; //x coordinate of 8 detector locations (cm)
    vector <float> ydet0; //y coordinate of 8 detector locations (cm)
    vector <float> pixx0; //x coordinate of 16 pixel locations (cm)
    vector <float> pixy0; //y coordinate of 16 pixel locations (cm)
    
    vector <float> pixxw; //width of pixel in x direction
    vector <float> pixyw; //width of pixel in y direction
    vector <vector <float> > areaDetector2D; //area of 16x16 pixels (in sq. cm) 
    vector <vector <float> > areaFullDetector2D; //area of 128x128 pixels (in sq. cm) 
    
    float areaGeom;       //Geometric area of all detectors
public:
    CztDetectorGeometry();
    /**
     * Defines the dimensiona and locations (lower left corners) of the detector modules
     * All lengths are in cm
     * @return 
     */
    int init_detectors();
    
    /**
     * Defines pixel coordinates of 16x16 pixels in a detector based on detector's [x,y] coordinate.
     * It defines their geometric width as well.
     * @param xdet: Detector X coordinate
     * @param ydet: Detector Y coordinate
     * @return 
     */
    int define_pixel_coordinates(float xdet, float ydet);
    

    
    /**
     * Generates 2D area map of all 128x128 pixels and stores them in areaFullDetector2D.
     * @return 
     */
    int get_area_map();
    
    /**
     * Displays detector coordinates for all 64 detectors and their complete geometric area
     * @return 
     */
    
    void display_detector_coordinates();
    
    void display_pixel_information(string infoType, float xdet=0.0, float ydet=0.0);
    
    
    //GETTER FUNCTIONS
    vector<float> get_xdet0(){return xdet0;}
    vector<float> get_ydet0(){return ydet0;}
    vector<float> get_pixx0(){return pixx0;}
    vector<float> get_pixy0(){return pixy0;}
    vector<float> get_pixxw(){return pixxw;}
    vector<float> get_pixyw(){return pixyw;}
    
    
    vector <vector <float> > get_area_full_detector(){ return areaFullDetector2D;}
    float get_geometric_area(){return areaGeom;}
    
};       
#endif	/* DETECTORGEOMETRY_H */
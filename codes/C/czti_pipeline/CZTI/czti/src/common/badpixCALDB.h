/* 
 * @file  badpixCALDB.h
 * @author Tanul Gupta
 * @date Created on May 18, 2015, 3:50 PM
 * @brief badpix CALDB file handler
 * @details This class contains function to read data and perform operation on CALDB badpix data.
 * @version 0.2.0
 */

#ifndef BADPIXCALDB_H
#define	BADPIXCALDB_H

//Include section
#include "glog/logging.h"
#include "utils.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fitsio.h>
#include <exception>
#include <map>
#include <algorithm>

// MACRO SECTION
struct BadpixTable{
    int quadID; //quadrant ID for this badpix table
    vector <unsigned char> detID;
    vector <unsigned char> pixID;
    vector <unsigned char> pixXquad;
    vector <unsigned char> pixYquad;
    vector <unsigned char> locX;
    vector <unsigned char> locY;
    vector <unsigned char> pixFlag;
    vector < vector <unsigned char> > vecBadpixMap; //badPixMap of size 64x64 to 
                                                    //be generated after reading badpix CALDB file.
    map<string, unsigned int> nBadpix; // 0-'GOOD'; 1-'SPECBAD'; 
                                       // 2-'FLICKERING'; 3-'NOISY'; 4-'DEAD' 
    
    
    int generate_badpix_info();
    /**
     * Creates 64x64 quadrant bad pixel image.
     * @return 
     */
    int create_quadrant_badPixMap();
    //GETTERS
    long get_nrows(){return detID.size();}
    //SETTERS
    void resize_table(long size);
    void pushback_pixFlag(unsigned char detid, unsigned char pixid, unsigned char pixxquad,
                        unsigned char pixyquad, unsigned char locx, unsigned char locy,
                        unsigned char pixflag);
    void set_quadID(int quadid){ this->quadID=quadid;}
    void reset();
    //Display
    void display();    
};

class Badpix {
private:
    // badpix[0] : for quadrant 0
    // badpix[1] : for quadrant 1
    // badpix[2] : for quadrant 2
    // badpix[3] : for quadrant 3
    // badpix[4] : for all quadrants
    vector <BadpixTable> badpix;
    vector <vector <unsigned char> > badpixMap; //badpix map for all 4 quadrants.

public:
    Badpix();
    /**
     * Reads CALDB badpixel file.
     * @param badpixFilename
     * @param extname
     * @return 
     */
    int read_badpix_file(string badpixFilename);
    
    /**
     * Writes level-2 bad pixel file
     * @param badpixFilename
     * @return 
     */
    int write_L2badpix_file(string badpixFilename);
    /**
     * Creates empty level-2 badpixel file
     * @param badpixFilename
     * @return 
     */
    int create_L2badpix_file(string badpixFilename, string badpixTemplate);
   
    int find_dead_noisy_pix(string dpiFilename, vector <vector <float> > emap);
    /**
     * Function to create and return 2-D bad pixel map of all quadrants.
     * @param status
     * @return 
     */
    void create_badPixMap();
    
    int display();
    
    //GETTERS
    /**
     * return badpix table for a particular quadrant number
     * @param quadNo
     * @return 
     */
    BadpixTable get_badpix_table(int quadNo){return badpix[quadNo];}
    
    /**
     * Return badpix table of all 4 quadrants
     * @return 
     */
    vector<BadpixTable> get_badpix_tables(){return badpix;}
    
    /**
     * Get badpix map for all quadrants
     * @return 
     */
    vector <vector <unsigned char> > get_badpix_map(){return badpixMap;}
    
    //SETTERS
    /**
     * @param bptables: badpix tables for all quadrants
     */
    int set_badpix_tables(vector <BadpixTable> bptables);
    
    /**
     * Clear badpix tables
     */
    void reset(){badpix.clear();}
    
};

/**
 * Compare multiple badpixel tables of a particular quadrant and creates a new badpix table
 * in which pixflag is the highest value for that pixel.
 * For Ex: if for a particular pixel, 3 badpixel tables are available; F1(flag in I table)=0,
 * F2=1, F3=4, then pixflag = Max(F1,F2,F3)=4.
 * CONSTRAINTS: all badpixel tables passed in the vector should have same quadrant id.
 * @param vecBadpixTable
 * @return 
 */
BadpixTable flag_badpix_quadrant(vector <BadpixTable> vecBadpixTable);

/**
 * Reads multiple badpixel files and evaluate highest pix flag value for every pixel (of
 * all 4 quadrants).
 * @param vecBadpixFiles
 * @return 
 */
vector<BadpixTable> flag_badpix_files(vector <string> vecBadpixFiles);


#endif	/* BADPIXCALDB_H */


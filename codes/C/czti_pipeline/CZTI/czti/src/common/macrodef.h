/* 
 * @file  macrodef.h
 * @author Tanul Gupta
 * @date Created on July 30, 2015, 2:51 PM
 * @brief Contains macros required in all pipeline modules.
 * @details Contains macros required in all pipeline modules.
 * @version 0.2.0
 */

 #ifndef MACRODEF_H
 #define MACRODEF_H

//****************Macro variable definitions****************
#define VERSION "2.1"
#define NAMESIZE 1024         //size of a filename/modulename

#define YES 1
#define NO 0
#define QUAD_MASKSIZE 8175
#define PIX_DIST_QUADRANTS 1850 //pixel distance between two quadrants at resolution of 0.02mm
#define MAX_KEYWORD_SIZE 1000
#define NO_QUAD_X 2
#define NO_QUAD_Y 2
#define NO_DET_X_PER_QUAD 4
#define NO_DET_Y_PER_QUAD 4
#define NO_DET_PER_QUAD 16 //number of detector per quadrant
#define NO_PIX_X_PER_DET 16
#define NO_PIX_Y_PER_DET 16
#define NO_PIX_ALL_QUADS 16384

#define NUMQUAD 4

#define XSIZE 128           //total pixels in x direction in CZTI
#define YSIZE 128
#define NUM_CHANNEL 1024      //NUmber of energy channels


#define XPIX_QUAD 64            //number of pixels in x axis in one quadrant
#define YPIX_QUAD 64

#define MAXCHANNEL 511      //Maximum channel number
#define MINCHANNEL 0       
#define EMAX   150          //maximum valid energy
#define EMIN   10          //minimum valid energy

#define ORIGIN "CZTI POC"
#define INSTRUME "CZTI"

#define EVENTDATA_HDU_BEGIN 2
#define EVENTDATA_HDU_END 5


//Pixel Quality 
#define GOODPIX 0
#define SPECBADPIX 1 //spectroscopically bad
#define FLICKPIX 2 //Flickering
#define NOISYPIX 3 //NOISY
#define DEADPIX 4 //DEAD/INACTIVE

#define BADPIX_EXT_CALDB "BADPIX"
#define BADPIX_VAL_KEYWORD "BADPIXL"  //in caldb file
#define GOODPIX_VAL_KEYWORD "GOODPIXL"
#define NOISYPIX_VAL_KEYWORD "NOISPIXL"
#define DEADPIX_VAL_KEYWORD "DEADPIXL"

#define ENV_PFILES "PFILES"

#define SHORT_STR_SIZE 25

#define NUM_HDU_EVENTFILE 7

//Macro constants for Detector
#define DET_MASKSIZE 1950

#define COMP_FIELD_SIZE 25

#define COMPCOLS 327
#define COMPROWS 8175

#define PIX_SIZ_MASK  0.002 //in cm
#define NPIXELS_GAP 125
#define CLOSE 0
#define PIXELS_PER_COL_DET 16
#define PIXELS_PER_ROW_DET 16

#define PIXELS_PER_ROW_QUAD 64
#define PIXELS_PER_COL_QUAD 64
#define NPIXELS_X 128
#define NPIXELS_Y 128

#define NDET_PER_ROW 8
#define NDET_PER_COL 8
#define NUM_DET_PER_QUAD 16     //Number of detectors per quadrant
#define NUM_DET_PER_QUAD_PER_ROW 4
#define NUM_DET_PER_QUAD_PER_COL 4

#define MASKHEIGHT 47.7         //in cm
#define DETECTORWIDTH 3.9       //in cm
#define TOTAL_DET_MODULES 64

#define AVG_PIXEL_WIDTH .24375  //in cm

#define XLDET 3.906 //x length in cm
#define YLDET 3.906 //y detector length in cm
#define XPITCH 4.15 //in cm
#define YPITCH 4.15 //in cm
#define XDMIN -18.203
#define YDMIN -18.203
#define ZMASK 47.7 //Height of mask (in cm)
#define CZTDET_THICKNESS 0.5 //CZT detector thickness (in cm)
#define PCB_THICKNESS 0.015 //PCB thickness (in cm)
#define TANMASK_THICKNESS 0.05 //Tantalum mass thickness (in cm)
#define BIG_PIXEL_WIDTH 0.246 // in cm
#define SMALL_PIXEL_WIDTH 0.231 // in cm

//mask geometry
#define XLMASK 3.9  //x length of mask per detector/module in cm
#define YLMASK 3.9  //y length of mask per detector/module in cm
#define XLMASK_PIXELS 1950 //in pixels
#define YLMASK_PIXELS 1950 //in pixels
#define XPITCH_MASK 4.15 //in cm
#define YPITCH_MASK 4.15 //in cm
#define XPITCH_MASK_PIXELS 2075 //in pixels
#define YPITCH_MASK_PIXELS 2075 //in pixels
#define XMIN_MASK -18.20 //in cm
#define YMIN_MASK -18.20 //in cm
#define NX_MASKPIXELS 18200 //number of sub-pixels in finer mask in x direction
#define NY_MASKPIXELS 18200 //number of sub-pixels in finer mask in y direction
#define SUB_PIXEL_WIDTH 0.002 //in cm
#define BIG_PIXEL_WIDTH_MASK_PIXELS 123 //in pixels
#define SMALL_PIXEL_WIDTH_MASK_PIXELS 114 //in pixels
#define DIFF_BIG_SMALL_PIXELS_IN_MASK 9 //in pixels
#define OPEN 1
#define CLOSE 0

//DELTA DEFINITIONS
#define DELTA_INSTRUMENT_TIME 0.1
//DEGREE CONVERTERS
#define TODEG 57.29578049
#define TORAD 0.017453292

//PHYSICAL CONSTANTS
#define RADEARTH 6378.135 //radius of earth (in km)

//TEMPLATE NAMES
#define EVTTEMPLATE "evtTemplate"
#define HDRTEMPLATE "hdrTemplate"
#define GTITEMPLATE "gtiTemplate"
#define BTITEMPLATE "btiTemplate"
#define BADPIXTEMPLATE "badpixTemplate"
#define SHADOWTEMPLATE "shadowTemplate"
#define CORMKFTEMPATE "correctedMKFTemplate"
#define IMAGETEMPLATE "imgTemplate2"
#define BUNCHTEMPLATE "bunchTemplate"
#endif /* MACRODEF_H */


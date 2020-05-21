/* 
 * File:   common.h
 * Author: preeti, Tanul Gupta
 *
 * Created on July 6, 2012, 5:31 PM
 * Modified on November 10, 2014 (by Tanul)
 * Modified on June 12, 2015 (by Tanul)
 * Modified on Sep 4, 2015 (by Tanul)
 */

#ifndef COMMON_H
#define	COMMON_H

//**********************Include files***********************
#include<pil.h>
#include "glog/logging.h"
#include "errorHandler.h"
#include<iostream>
#include<fitsio.h>
#include<cstdlib>
#include<cstring>
#include<string>
#include<stdarg.h>
#include<cmath>
#include<dirent.h>
#include<unistd.h>
#include<strings.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<vector>
#include<map>
#include<numeric>
#include<algorithm>
#include<cstdio>
#include<iomanip>
#include<sstream>
#include "macro_definitions.h"
#include <string>
#include "cztstring.h"
#include "macrodef.h"
#include "errorHandler.h"

using namespace std;


template <class T>
void readFitsKey(fitsfile *fptr, int DTYPE, string keyname, T *value,
        string *comment){
    ErrorHandler errHandler;
    int status=0;
    char chrComment[PIL_LINESIZE];
    fits_read_key(fptr, DTYPE, (char*) keyname.c_str(), value, chrComment, 
            & status);

    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading fits key " + keyname;
        throw errHandler;
    }
}

template <class T>
void readFitsKeyStr(fitsfile *fptr, string keyname, T *value,
        string *comment){
    ErrorHandler errHandler;
    int status=0;
    char chrComment[PIL_LINESIZE];
    char keyValue[PIL_LINESIZE];
    fits_read_key(fptr, TSTRING, (char*) keyname.c_str(), keyValue, chrComment, 
            & status);
    
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading fits key.";
        //throw errHandler;
    } else {
        *value = (string) keyValue;
    }
}

int copyKeywords(fitsfile *fptr1,fitsfile *fptr2,int n,...);
//------------------------------------------------------------------

/**
 * Function to copy user keywords from input fitsfile to output fitsfile.
 * It works on the current HDUs of both fits file 
 * @param fin
 * @param fout
 * @return 
 */
int copyUserKeyWords(fitsfile *fin,fitsfile *fout);

/**
 * Function to copy a number of keywords from one fits file extension to another fits file extension.
 * @param fin
 * @param fout
 * @param inHDU
 * @param outHDU
 * @param vecKeywords
 * @return 
 */
int copyUserKeyWords(fitsfile *fin, fitsfile *fout, 
        string inHDU, string outHDU, vector <string> vecKeywords);
/**
 * Function to print error based on fits error code
 * @param status
 * @param errstring
 */
//void printerror(int status,string errstring="Error");
//------------------------------------------------------------------

/**
 * Function to join two or more strings.
 * @param outstring : the output string. It should have enough space to
 * hold the concatenated string
 * @param n : number of strings to be joined
 * @param ... : strings to join
 * @return  new concatenated string
 */
void joinStrings(char *outstring,int n,...);
//-------------------------------------------------------------------

/**
 * Function to delete non empty directory recursively
 * @param dir
 * @return 
 */
int deleteDir(char *dir);
//-------------------------------------------------------------------

/**
 * Function to filter out . and .. directories from output of scandir function
 * @param dptr
 * @return returns 0 for . and .. directories and 1 for others
 */
int defaultfilter(const struct dirent *dptr);
//-------------------------------------------------------------------

/**
 * Parses the string depending on the delimiter
 * @param str
 * @param delim
 * @param substr
 * @return 
 */
int parseString(string str,char delim,vector<string> &substr);
//-------------------------------------------------------------------

//-------------------------------------------------------------------

/**
 * Function to get start and end of bins from an input string
 * Ex. - if ebins supllied is like - "1-4,5-8,9-12" then getbins will store them in two arrays, namely start and end in the
 * following manner- start=[1,5,9] & end=[4,8,12]
 * @param ebins
 * @param start
 * @param end
 * @param nebins
 * @return
 */
int getbins(char *ebins,double *start,double *end,int nebins);
//-------------------------------------------------------------------

template <class T>
int getVecBins(string bins, vector<T> *startVec, vector<T> *stopVec, int nbins) {
    int status=0;
    int ibins=0;
    int pos=0;
    string temp, temp2;
    if (bins=="-"){
        return EXIT_SUCCESS;
    } else {
        (*startVec).resize(nbins, 0.0);
        (*stopVec).resize(nbins, 0.0);
        for(ibins=0; ibins<nbins; ibins++) {
            pos = bins.find('-', 0);
            temp = bins.substr(0, pos);
            temp2 = bins.substr(pos + 1);
            bins.clear();
            bins.assign(temp2);
            (*startVec)[ibins] = (T) atof(temp.c_str());
            temp.clear();
            temp2.clear();

            pos = bins.find(',', pos);
            temp = bins.substr(0, pos);
            temp2 = bins.substr(pos + 1);
            bins.clear();
            bins.assign(temp2);
            (*stopVec)[ibins] = (T) atof(temp.c_str());
            temp.clear();
            temp2.clear();
        }
    }
    return EXIT_SUCCESS;
}

/**
 * Function to get number of bins separated by comma from input string or as 
 * first line in an ascii input file  
 * @param bins
 * @return 
 */
int get_nbins(char *bins);
//-------------------------------------------------------------------

/**
 * Function to dynamically allocate memory to 2D array 
 * @param h
 * @param w
 * @return 
 */
template<class T> T **allocateMemory(long h,long w);
//-------------------------------------------------------------------

/**
 * Function to initialize each element of array with val 
 * @param array
 * @param size
 * @param value
 * @return 
 */
template<class T> int initArray(T *array,long size, T value);
//-------------------------------------------------------------------

/**
 * Function to free dynamically allocated 2D memory
 * @param array
 * @param height
 * @param width
 */
template<class T> void freeMemory(T **array,long height);
//-------------------------------------------------------------------

/**
 * Function to compute rms of an array
 * @param arr
 * @param size
 * @return 
 */
template<class T> double getrms(T *arr,long size);
//-------------------------------------------------------------------

/**
 * Function to compute mean of an array
 * @param arr
 * @param size
 * @return 
 */
template<class T> double getmean(T *arr,long size);
//-------------------------------------------------------------------

/**
 * Function to compute maximum value from an array
 * @param array
 * @param size
 * @param maxval
 * @return 
 */
template<class T> long max(T *array,long size,T *maxval);
//-------------------------------------------------------------------

/**
 * Function to compute minimum value from an array
 * @param array
 * @param size
 * @param minval
 * @return 
 */
template<class T> long min(T *array,long size,T *minval);
//-------------------------------------------------------------------

//****************Function Definitions************************
template<class T>
T **allocateMemory(long h,long w){
    T **array;
    array=new T*[h];
    if(array==NULL)  return NULL;
    for(long i=0;i<h;i++){
        array[i]=new T[w]; 
        if(array[i]==NULL)
            return NULL;
    }
    return array;
}
//-------------------------------------------------------------------
template<class T>
int initArray(T *array,long size, T value){
    for(long i=0;i<size;i++){
        array[i]=value;
    }
    return (EXIT_SUCCESS);
}
//-------------------------------------------------------------------
template<class T>
void freeMemory(T **array,long height){
    for(long i=0;i<height;i++){
        delete[] array[i];
    }
    delete[] array;
}
//-------------------------------------------------------------------
template<class T>
double getrms(T *arr,long size){
    long i;
    double sum=0,rms=0;
    for(i=0;i<size;i++){
        sum=sum+arr[i]*arr[i];
    }
    //sum=sqrt(sum);
    rms=(double)sqrt(sum/(double)size);
    return rms;
}
//-------------------------------------------------------------------             
template<class T>
double getmean(T *arr,long size){
    long i;
    double sum=0,mean=0;
    for(i=0;i<size;i++){
        sum=sum+arr[i];
    }
    mean=sum/(double)size;
    return mean;
}
//-------------------------------------------------------------------
template<class T>
long max(T *array,long size,T *maxval){
    long index=0;
    if(size>0){
        *maxval=array[0];
        for(int i=0;i<size;i++){
            if(*maxval<array[i]){
                *maxval=array[i];
                index=i;
            }
        }
        return index;
     }
    *maxval=0;
    return (EXIT_FAILURE);
}
//-------------------------------------------------------------------
template<class T>
long min(T *array,long size,T *minval){
    long index=0;
    if(size>0){
        *minval=array[0];
        for(int i=0;i<size;i++){
            if(*minval>array[i]){
                *minval=array[i];
                index=i;
            }
        }
        return index;
     }
    *minval=0;
    return (EXIT_FAILURE);
}

//---------Quaternion related structure and functions---------
struct Q{
    double q1,q2,q3,q4;   //quaternion representation  q1-scaler
    double theta,x,y,z;   //for axis angle representation
    double mod;             //modulus of quaternion
    bool norm,AxisAngle,qflag;           
    Q(); 
    Q(double a, double b, double c, double d);
    Q(const Q &q);
    void readQ();
    void readAxisAngle();
    int getAxisAngle();
    int getQuat();
    void display();
    void normalize();   //to make modulus 1
 };
//-------------------------------------------------------------------
struct Axis{
    double x,y,z,mod;
    bool norm;
    Axis(){ x=0; y=0; z=0; norm=false;}
    Axis(double a, double b, double c){
        x=a; y=b; z=c; norm=false;
    }
    double getMod();
    void normalize();   //to make modulus one
    void display();
};
//-------------------------------------------------------------------

/**
 * Function to get the product of two float matrices A and B. The product is given
 * in matric C. Size of A if axb, size of B is bxc and size of C is axc
 * @param A
 * @param B
 * @param C
 * @param a
 * @param b
 * @param c
 */
void matrix_product(float **A,float **B,float **C,int a,int b,int c);
//-------------------------------------------------------------------

/**
 * Function to multiply two quaternions
 * Q1 x Q2 = Q3
 * @param Q1
 * @param Q2
 * @param Q3
 */
void quaternion_product(Q &Q1, Q &Q2, Q &Q3);
//-------------------------------------------------------------------

/**
 * Function to return inverse of a quaternion
 * @param q
 * @return 
 */
void Inverse(Q &q,Q &inverseq);
//-------------------------------------------------------------------

/**
 * Transform the axis to a different frame using the given quaternion
 * @param axis
 * @param q
 * @param axisrot
 * @return 
 */
int rotate(Axis &axis, Q &q, Axis &axisrot);
//-------------------------------------------------------------------

/**
 * Compute Camera coordinates (theta_x,theta_y) using transformation qc from 
 * average aspect file and RA DEC
 * @param aspectfile : Average aspect file
 * @param RA : RA in radians
 * @param DEC : Declination in radians
 * @param theta_x : In radians
 * @param theta_y : In radians
 * @return  
 */
int getCC(char *aspectfile,double RA,double DEC,double *theta_x,double *theta_y);
//-------------------------------------------------------------------------------

/**
 * Computes RA & DEC (in radians) using qc from average aspect file and Camera 
 * Coordinates.
 * @param aspectfile: Average aspect file
 * @param thetax: camera coordinate x in radians
 * @param thetay: camera coordinate y in radians
 * @param RA: RA in radians
 * @param DEC: DEC in radians
 */
int getRaDec(char *aspectfile,double thetax,double thetay,double *RA, double *DEC);

/**
 * Function to check the validity of range of values for the given column in a given hdu
 * @param fptr
 * @param hdunum
 * @param colname
 * @param val1
 * @param val2
 * @return 0 for valid, 1 for invalid, -1 for error
 */
template<class T>
int checkRangeValidity(fitsfile *fptr,int hdunum,char *colname,T *val1,T *val2,int num, bool *valid){
    int status=0;
    int hdutype;
    for(int i=0;i<num;i++){
        if(val1[i]>=val2[i]){
            cout<<"\n***Improper range given***\n";
            return (EXIT_FAILURE);
        }
    }
    //cout<<"\nInput Range:"<<val1<<"-"<<val2;
    fits_movabs_hdu(fptr,hdunum,NULL,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    fits_get_hdu_type(fptr,&hdutype,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    if(hdutype==IMAGE_HDU){
        cout<<"\n***Function valid only for Table HDU-checkRangeValidity()***\n";
        return (EXIT_FAILURE);
    }
  
    int colnum;
    fits_get_colnum(fptr,CASEINSEN,colname,&colnum,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    cout<<"\nColumn number for "<<colname<<" is "<<colnum;
    int typecode=0;
    fits_get_eqcoltype(fptr,colnum,&typecode,NULL,NULL,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    
    long nrows;
    fits_get_num_rows(fptr,&nrows,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    
    switch(typecode){
        case TBIT:   
        case TBYTE:  
        case TSBYTE:  
        case TUSHORT:
        case TSHORT:  
        case TFLOAT:  
        case TINT:
        case TUINT:   
        case TDOUBLE:
        case TULONG:
        case TLONGLONG:
        case TLONG:    
            break;
        default: cout<<"\nThis function does not supports typecode "<<typecode;
            return (EXIT_FAILURE);
    }
        
    double *array=(double *)malloc(nrows*sizeof(double));
    if(array==NULL) {
        cout<<"***Out of memory error-checkRangeValidity()***\n";
        return (EXIT_FAILURE);
    }
    
    if(nrows>0){
        long firstrow=1;
        long firstelem=1;
        cout<<"\nReading column number "<<colnum;
        fits_read_col(fptr,TDOUBLE,colnum,firstrow,firstelem,nrows,NULL,array,NULL,&status);
        if(status){
            fits_report_error(stderr,status);
            return (EXIT_FAILURE);
        }
        
        T maxval=0,minval=0;
        max(array,nrows,&maxval);
        min(array,nrows,&minval);
        free(array);
        for(int i=0;i<num;i++){
            if(val1[i]>maxval || val2[i]<minval)
                valid[i]=false;
            
            else
                valid[i]=true;
        }
    }
    return (EXIT_SUCCESS);
               //invalid as there is no data
 }
//-------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
/**
 * Function to write history to all HDUs of the file
 * @param filename
 * @param vhistory
 * @return 
 */
int writeHistory(char *filename,vector<string> &vhistory);
//---------------------------------------------------------------------------------------

/**
* Template function to evaluate min and maximum value of a vector of any data type.
* @param first : index of first value in vector
* @param last : index of last value in vector
* @param min : stores the minimum value found in the index range specified
* @param max : stores the maximum value found in the index range
*/
template<class Iter_T, class Value_T> 
void compute_min_and_max(Iter_T first, Iter_T last, Value_T &min, Value_T &max) {
    min = *min_element(first, last);
    max = *max_element(first, last);
}
//for date , checksum, origin, creator

//Reading Header keywords
template <class T>
int readKey(fitsfile *fptr, int DATATYPE, string keyname, T *value){
    ErrorHandler errHandler;
    int status =0;
    fits_read_key(fptr, DATATYPE, keyname.c_str(), value, NULL, &status);
    if(status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading key " + keyname;
        throw errHandler;
    }
    
    return status;
}

//Update Header key

template <class T>
int updateKey(fitsfile *fptr, int DATATYPE, string keyname, T *value, string comment="") {
    ErrorHandler errHandler;
    int status = 0;
    fits_update_key(fptr, DATATYPE, keyname.c_str(), value, (char*) comment.c_str(), &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in updating key " + keyname;
        throw errHandler;
    }

    return status;
}


/**
 * Function to add keywords DATE, CREATOR, ORIGIN and CHECKSUM to all hdus of file
 * @param filename
 * @param creator
 * @return 
 */
int updateKeywords(char *filename,char *creator);

int updateHdrTime(string filename, string extname, string colname);

   
/**
 * Update/Write key and corresponding values
 * The function expects to be run multiple times for different data types.
 * @param filename
 * @param extname
 * @param keys
 * @param values
 * @param DATA_TYPE data type of key values to be written (ex. TFLOAT, TINT)
 * @return 
 */
template<class T>
int updateKeys(string filename, string extname, vector<string> keys, vector<T> values, int DATA_TYPE){
    int status=0;
    int ikey=0, nkeys=0;
    fitsfile *fptr;
    char keyValue[PIL_LINESIZE];

    fits_open_file(&fptr, (char*) filename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeys()";
        fits_report_error(stderr, status);
        return status;
    }

    fits_movnam_hdu(fptr, BINARY_TBL, (char*) (extname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << extname << " of " <<
                filename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    nkeys=values.size();

    for (ikey = 0; ikey < nkeys; ikey++) {
        fits_update_key(fptr, DATA_TYPE, (char*) (keys[ikey]).c_str(), &values[ikey], NULL, &status);
        if (status) {
            LOG(WARNING) << "Unable to update key " << keys[ikey];
        }
    }

    
    fits_close_file(fptr, &status);
    if(status){
        LOG(ERROR) << "Error in closing file " << filename;
        return EXIT_FAILURE;
    }
    
    return status;
}
/**
 * Update/Write key and corresponding values
 * The function expects to be run multiple times for different data types.
 * @param fptr
 * @param extname
 * @param keys
 * @param values
 * @param DATA_TYPE data type of key values to be written (ex. TFLOAT, TINT)
 * @return 
 */
template<class T>
int updateKeys(fitsfile *fptr, string extname, vector<string> keys, vector<T> values, int DATA_TYPE){
    int status=0;
    int ikey=0, nkeys=0;
    char filename[FLEN_FILENAME];
    char keyValue[PIL_LINESIZE];

    fits_file_name(fptr, filename, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeys()";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) (extname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << extname << " of " <<
                filename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    nkeys = keys.size();
    LOG(INFO) << "Updating keys in file: " << filename;
    for (ikey = 0; ikey < nkeys; ikey++) {
        fits_update_key(fptr, DATA_TYPE, (char*) (keys[ikey]).c_str(), &values[ikey], NULL, &status);
        if (status) {
            LOG(WARNING) << "Unable to update key " << keys[ikey];
        }
    }

    
    return status;
}
/**
 * Update/Write key and corresponding values
 * The function expects to be run multiple times for different data types.
 * @param filename
 * @param extname
 * @param keys
 * @param values
 * @param DATA_TYPE data type of key values to be written (ex. TFLOAT, TINT)
 * @return 
 */
template<class T>
int updateKeys(string filename, int extnum, vector<string> keys, vector<T> values, int DATA_TYPE){
    int status=0;
    int ikey=0, nkeys=0;

    fitsfile *fptr;

    fits_open_file(&fptr, (char*) filename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeys()";
        fits_report_error(stderr, status);
        return status;
    }

    fits_movabs_hdu(fptr, extnum, NULL, &status);
    if(status) {
        LOG(ERROR) << "Error in moving to HDU " << extnum << " of " <<
                filename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    nkeys=values.size();

        for (ikey = 0; ikey < nkeys; ikey++) {
            fits_update_key(fptr, DATA_TYPE, (char*) (keys[ikey]).c_str(), &values[ikey], NULL, &status);
            if (status) {
                LOG(WARNING) << "Unable to update key " << keys[ikey];
            }
        }
    
    fits_close_file(fptr, &status);
    if(status){
        LOG(ERROR) << "Error in closing file " << filename;
        return EXIT_FAILURE;
    }
    
    return status;
}

template <class T>

//-----------------------------------------------------------------------------------------

/**
 * Function to create full image from data of four different quadrants where the input quadrant
 * arrays are each of single dimension 64x64 and the output array is a single dimension array of 128x128 size.
 * @param q1
 * @param q2
 * @param q3
 * @param q4
 * @param fullImage
 * @return 
 */
/****************************************************************************************************

 ****************************************************************************************************/

int makeFullImage(long q0[XPIX_QUAD*YPIX_QUAD],long q1[XPIX_QUAD*YPIX_QUAD],
        long q2[XPIX_QUAD*YPIX_QUAD],long q3[XPIX_QUAD*YPIX_QUAD],long fullImage[XSIZE*YSIZE]);
//---------------------------------------------------------------------------------------------

template<class T>
int makeFullImage(T *q0,T *q1,T *q2,T *q3,T *fullImage){
     /**
     * Starting co-ordinates for quadrants are-
     * Q0   (0,0)     (0,0)
     * Q1   (64,0)    (XPIX_QUAD,0)
     * Q2   (64,64)   (XPIX_QUAD,YPIX_QUAD) 
     * Q3   (0,64)    (0,YPIX_QUAD)
     */
    //for 1st quadrant
    //loop cariables k,l for fullImage
    //loop variables i,j for qs
    for(int i=0;i<XSIZE*YSIZE;i++)
        fullImage[i]=0;
    
    for(int i=0,k=0;i<YPIX_QUAD;i++,k++){
        for(int j=0,l=0;j<XPIX_QUAD;j++,l++){
            fullImage[k*XSIZE+l]=q0[i*XPIX_QUAD+j];
        }
    }
    //for 2nd quadrant
    for(int i=0,k=0;i<YPIX_QUAD;i++,k++){
        for(int j=0,l=XPIX_QUAD;j<XPIX_QUAD;j++,l++){
            fullImage[k*XSIZE+l]=q1[i*XPIX_QUAD+j];
        }
    }
    //for 3rd quadrant
    for(int i=0,k=YPIX_QUAD;i<YPIX_QUAD;i++,k++){
        for(int j=0,l=XPIX_QUAD;j<XPIX_QUAD;j++,l++){
            fullImage[k*XSIZE+l]=q2[i*XPIX_QUAD+j];
        }
    }
    //for 4th quadrant
   for(int i=0,k=YPIX_QUAD;i<YPIX_QUAD;i++,k++){
        for(int j=0,l=0;j<XPIX_QUAD;j++,l++){
            fullImage[k*XSIZE+l]=q3[i*XPIX_QUAD+j];
        }
    }  
    return (EXIT_SUCCESS);
}
//----------------------------------------------------------------------------------------------
/**
 * Function to merge two or more image counts. Size gives the total number of pixels 
 * in the image. Image argument  must be given as type long *
 * @param long*
 * @param size
 * @param n
 * @param ...
 * @return 
 */
int mergeImages(long *out,int size,int n,...);




//functions that maps quadID, detID, pixID to get the image coordinates in x,y taking all quadrants
int getxy(int quadID,int detID,int pixID, unsigned char *x,unsigned char *y);

int getRaDec(char *file,int hdunum,double *RA,double *DEC);

int writeImg(char *file, long *img,int m,int n);  //just to check;
int writeImg(char *file, float *img,int m,int n);  //just to check
int writeImgD(char *file, double *img, int m, int n); //just to check;
int writeArray(char *file,float **array,int m,int n);


/**
 * Function to return quadrant id in quad and x,y location in that quadrant given
 * x,y in full 128 x 128 array
 * @param x
 * @param y
 * @param qx
 * @param qy
 * @param quad
 * @return 
 */
int getqxqy(int x,int y,int *qx,int *qy,int *quad);

/**
 * Function to return x,y in full 128 x 128 array when x,y is given wrt one quadrant 
 * @param qid
 * @param qx
 * @param qy
 * @param x
 * @param y
 * @return 
 */
int getxy4mqxqy(int qid,int qx,int qy,int *x,int *y);
 
void checkPFILESenv();

void checkParFile(char *modulename);

/**
 * Function to read Quantum efficiency values from Gain Offset QE file to the array qearray of size XSOZE*YSIZE
 * @param qefile
 * @param qearray
 * @return 
 */
int readQE(string qefile,float *qearray);



/******************************************************************************************************
 *Function to convert quadrant number values into its extension name in event file.
 * @param quad_no - quadrant number input to the function.
 * @return HDUname - stores the extension name for specified quadrant number. 
 * E.x. - 0 - 'q0', 1 - 'q1', 2 - 'q2, 3 - 'q3'
 ******************************************************************************************************/
void quadToHDU(int quad_no, char* HDUname);

/**
 * Converts HDU name into corresponding quadrant number.
 * @param HDUname
 * @param quadNo
 */
void HDUtoQuad(string HDUname, unsigned char &quadNo);
/******************************************************************************************************
 *Function to convert quadrant number values into its extension name in dpi file.
 * @param quad_no - quadrant number input to the function.
 * @return HDUname - stores the extension name for specified quadrant number. 
 * E.x. - 0 - 'q0', 1 - 'q1', 2 - 'q2, 3 - 'q3'
 ******************************************************************************************************/
void DPIquadToHDU(int quad_no, char* HDUname);




int breakDPIarray(float *dpi, float **quad_dpi, int full_XSIZE=128, int full_YSIZE=128, int quad_XSIZE=64, int quad_YSIZE=64);



//>>>>>>>>>>>>>>>>>> CONVERSIONS
template<class T> double angleConverter(T val_angle, char* input_unit, char* output_unit);

template<class T> double angleConverter(T val_angle_in, char* input_unit, char* output_unit){
    double val_angle_out=0;
    if(strcmp(input_unit, "arcmin")==0){
        if(strcmp(output_unit, "radians")==0){
            val_angle_out = ((double)val_angle_in/60)*M_PI/180;
//            cout << val_angle_out; //DEBUG
            return val_angle_out;
        }
    }
    return val_angle_out;
}

//>>>>>>>>>>>>>>>>>>> CONVERSIONS END

/**
 * This function is general function to handle while performing cfitsio functions
 * on fits data.
 * @param status: if status=0 then it does nothing else it generates an error with
 *                the output provided by user in error_string.
 * @param error_string: string explaining the error (provided by developer)
 * @return 
 */
int report_error(int status, string error_string);

/**
 * This function returns extension name of quadrant in event file
 * @param quadNo: quadrant number
 * @return 
 */
string quadExtName(int quadNo); 

/**
 * 
 * @param x: input X
 * @param x1
 * @param x2
 * @param y1
 * @param y2
 * @param y: outputY
 * @return 
 */
//int find_interpolated_value(float x, float x1, float x2, float y1, float y2, float &y);
/**
 * Function to interpolate and get value by linearly fitting two input values.
 *     y2-y1 
 * y = ----- * x2-x1  + y1
 *     x2-x1     
 * @param x1
 * @param y1
 * @param x2
 * @param y2
 * @param x
 * @param y
 * @return 
 */
template <class T> int find_interpolated_value(T x, T x1, T x2, T y1, T y2, T &y){
    T slope;
    if (x2==x1){
        LOG(ERROR) << "No interpolation possible as x2=x1 i.e. line is parallel to y axis.";
        return EXIT_FAILURE;
    }
    slope = (y2-y1)/(x2-x1);
    y = (slope*(x-x1)) + y1;
    
    return EXIT_SUCCESS;
}

/**
 * Checks whether a particular pixel is present in quads (supplied by user);
 * returns true if pixel is present otherwise false
 * @param locx
 * @param locy
 * @param quadsToProcess
 * @param quadID- Quadrant ID in which pixel is present.
 * @param pixx- location x of pixel in a single quadrant
 * @param pixy- location y of pixel in a single quadrant
 * @return 
 */
bool is_pix_present_in_quads(unsigned char locx, unsigned char locy, vector<int> quadsToProcess, int &quadID, int &pixx, int &pixy);
/**
 * Generates detX, detY from detID, pixID & quadID
 * @param detID: detector ID
 * @param pixID: pixel ID
 * @param quadid: quadrant ID
 * @param detx: detector X coordinate
 * @param dety: detector Y coordinate
 * @return 
 */
template <class T> int generate_detx_dety(T detID, T pixID, T quadid, T &detx, T &dety) {
    //T x, y, px, py;
//    x = (detID % 4)*16;
//    y = (3 - (detID / 4))*16;
//    px = pixID % 16;
//    py = 16 - (ceil((double) (pixID + 1) / 16.0));
//    if (quadid == 1 || quadid == 2) {
//        detx = 63 - (x + px);
//        dety = 63 - (y + py);
//    } else {
//        detx = x + px;
//        dety = y + py;
//    }

    T px, py, dx, dy;
    dx = detID % 4;
    dy = detID / 4;
    px = pixID % 16;
    py = pixID / 16;
    detx = (dx * 16 + px);
    dety = dy * 16 + py;
    if(quadid==0 ||quadid==3){
        dety=63-dety;
    } else {
        detx=63-detx;
    }
    
    return EXIT_SUCCESS;
}

template <class T> int generate_locx_locy(T detID, T pixID, T quadid, T &detx, T &dety, T &locx, T &locy){
    T xoff=0, yoff=0;
    if(generate_detx_dety(detID, pixID, quadid, detx, dety)){
        LOG(ERROR) << "Error in generating actual pixel coordinates for quadrant " << quadid;
        return EXIT_FAILURE;
    }
    switch (quadid) {
        case 0:
            xoff = 0;
            yoff = 64;
            break;
        case 1:
            xoff = 64;
            yoff = 64;
            break;
        case 2:
            xoff = 64;
            yoff = 0;
            break;
        case 3:
            xoff = 0;
            yoff = 0;
            break;
    }
    locx = detx + xoff;
    locy = dety + yoff;
    
    return EXIT_SUCCESS;
}
template <class T> int generate_locx_locy(T detx, T quadid, T dety, T &locx, T &locy){
    T xoff=0, yoff=0;
    switch (quadid) {
        case 0:
            xoff = 0;
            yoff = 64;
            break;
        case 1:
            xoff = 64;
            yoff = 64;
            break;
        case 2:
            xoff = 64;
            yoff = 0;
            break;
        case 3:
            xoff = 0;
            yoff = 0;
            break;
    }
    locx = detx + xoff;
    locy = dety + yoff;
    
    return EXIT_SUCCESS;
}

/*
template <class T> int generate_detx_dety_old(T detID, T pixID, T quadid, T &detx, T &dety) {
    T px, py, dx, dy;
    dx = detID % 4;
    dy = detID / 4;
    px = pixID % 16;
    py = pixID / 16;
    detx = dx * 16 + px;
    dety = dy * 16 + py;
    
    return EXIT_SUCCESS;

}
*/

/**
 * Creates empty fits file from template file
 * @param outFilename
 * @param outTemplate
 * @return 
 */
int create_empty_fitsfile(string outFilename, string outTemplate);
/**
 * This function generates full path of Template file by reading the value of environment
 * variable CZTI_templates
 * @param filename: name of template file
 * @return complete path to template file
 */
string template_full_path_generator(string filename);


template<class T> int linearize_2Dvector(vector <vector <T> > vector_2D, vector<T> *vector_1D){
    long ix, iy =0;
    long nx, ny =0;
    ErrorHandler errHandler;
    ny =vector_2D.size();
    if(ny==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Empty input 2D vector.";
        throw errHandler;
    }
    nx = vector_2D[1].size();
    (*vector_1D).resize(ny*nx, 0.0);
    for(iy=0; iy<ny; iy++){
        for(ix=0; ix<nx; ix++){
            (*vector_1D)[iy*nx + ix] = vector_2D[iy][ix];
        }
    }
    
    return EXIT_SUCCESS;
}

/**
 * Function to create a 2D vector from a 1-Dimensional vector
 * @param y: Vector containing row numbers
 * @param x: Vector containing column numbers
 * @param input_vector: vector to be rearranged
 * @param rearranged_vector: rearranged vector
 * @return 
 */
template <class Index_T, class Array_T>
int create_2d_vector(vector<Index_T> y, vector<Index_T> x, vector<Array_T> input_vector, vector <vector <Array_T> > *rearranged_vector){
    int i,j=0;
    int status=0;
    for (i=0; i<y.size(); i++){
        (*rearranged_vector)[y[i]][x[i]] = input_vector[i];   
    }
    return status;
}


/*template <class Index_T, class Array_T>
int create_2d_vector_to_array(vector<Index_T> y, vector<Index_T> x, vector<Array_T> input_vector,Array_T *(**rearranged_array)){
    int i,j=0;
    int status=0;
	Array_T *input_vector_to_array;
    for (i=0; i<y.size(); i++){
        input_vector_to_array[i]= input_vector[i]; 
	    (*rearranged_array)[y[i]][x[i]] = input_vector_to_array[i]; 	
    }
    return status;
}*/


/**
 * Rearrange Quad to properly place it in a full image based on it's quadrant number.
 * @param fullImage
 * @param quadImage
 * @param quadNo
 * @return 
 */
template <class T> int rearrange_quads(vector <vector <T> > *fullImage, vector <vector <T> > quadImage, int quadNo){
    int ix, iy=0;
    int xoff, yoff=0;
    int status=0;
    
    switch (quadNo) {
        case 0:
            xoff=0;
            yoff=64;
            break;
        case 1:
            xoff=64;
            yoff=64;
            break;
        case 2:
            xoff=64;
            yoff=0;
            break;
        case 3:
            xoff=0;
            yoff=0;
            break;
    }

    for (iy = 0; iy < YPIX_QUAD; iy++) {
        for (ix = 0; ix < XPIX_QUAD; ix++) {
            (*fullImage)[iy+yoff][ix+xoff] = quadImage[iy][ix];
        }
    }    
    return status;
}

/**
 * Rearrange finer mask quads to properly place it in a full image based on it's quadrant number.
 * @param fullImage
 * @param quadImage
 * @param quadNo
 * @return 
 */
template <class T> 
int rearrange_finer_mask_quads(vector <vector <T> > *fullImage, vector <vector <T> > quadImage, int quadNo) {
    int ix, iy = 0;
    int xoff, yoff = 0;
    int status = 0;

    switch (quadNo) {
        case 0:
            xoff = 0;
            yoff = 8175 + PIX_DIST_QUADRANTS;
            break;
        case 1:
            xoff = 8175 + PIX_DIST_QUADRANTS;
            yoff = 8175 + PIX_DIST_QUADRANTS;
            break;
        case 2:
            xoff = 8175 + PIX_DIST_QUADRANTS;
            yoff = 0;
            break;
        case 3:
            xoff = 0;
            yoff = 0;
            break;
    }

    for (iy = 0; iy < QUAD_MASKSIZE; iy++) {
        for (ix = 0; ix < QUAD_MASKSIZE; ix++) {
            (*fullImage)[iy + yoff][ix + xoff] = quadImage[iy][ix];
        }
    }
    return status;
}

template <class T>
int rearrange_detectors(vector <vector <T> > *quadImage, vector <vector <T> > detImage, int detNo ){
    int ix, iy=0;
    int xoff, yoff =0;
    int status=0;
    int row=0, col=0;
    
    row=detNo/4;
    col=detNo%4;
    
    xoff = col*16;
    yoff = row*16;
    
    for(iy=0; iy<NO_PIX_Y_PER_DET; iy++){
        for(ix=0; ix<NO_PIX_X_PER_DET; ix++){
            (*quadImage)[iy+yoff][ix+xoff] = detImage[iy][ix];
        }
    }
    
    return status;
}

template <class T>
string itoa(T number, int precision=5){
    ostringstream ss;
    ss << setprecision(precision) << number;
    return ss.str();
}


//DEBUG FUNCTIONS
template <class T> int write_img(vector <vector <T> > image, string outFileName, int IMAGE_TYPE, int DATA_TYPE){
    fitsfile *fimg;
    vector <T> image1D;
    int status=0;
    int bitpix = IMAGE_TYPE;
    int naxis=2;
    long naxes[2];
    long fpixel[2];
    long imgSize;
    fpixel[0] = fpixel[1]=1;
    naxes[1] = image.size();
    naxes[0] = 0;
    if(image.size()){
    naxes[0] = image[0].size();
    }
    imgSize=naxes[0]*naxes[1];

    DLOG(INFO) << "[DEBUG] Writing image with dimensions " << naxes[0] << " x " << naxes[1];
    if(imgSize==0){
        LOG(ERROR) << "Image data empty.";
        return EXIT_FAILURE;
    }
    fits_create_file(&fimg, (char*) outFileName.c_str(), &status);
    if (status) {
        LOG(ERROR) << "Error in creating image file " << outFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_create_img(fimg, bitpix, naxis, naxes, &status);
    if (status) {
        LOG(ERROR) << "Error in creating new image HDU to hold image data";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    DLOG(INFO) << "Linearizing 2D image for writing in fits";
    
    linearize_2Dvector(image, &image1D);
    fits_write_pix(fimg, DATA_TYPE, fpixel, imgSize, image1D.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing image data";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fimg, &status);
    if (status) {
        LOG(ERROR) << "Error in closing output image file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}

/**
 * Write image extension
 * @param image
 * @param outFileName
 * @param IMAGE_TYPE
 * @param DATA_TYPE
 * @return 
 */
template <class T> int write_img(vector <vector <T> > image, fitsfile *fimg, string extname, int IMAGE_TYPE, int DATA_TYPE){
    vector <T> image1D;
    int status=0;
    string errorMsg="";
    int bitpix = IMAGE_TYPE;
    int naxis=2;
    long naxes[2];
    long fpixel[2];
    long imgSize;
    fpixel[0] = fpixel[1]=1;
    naxes[1] = image.size();
    naxes[0] = 0;
    if (image.size()) {
        naxes[0] = image[0].size();
    }
    imgSize=naxes[0]*naxes[1];

    if (imgSize == 0) {
        LOG(ERROR) << "Image data empty.";
        return EXIT_FAILURE;
    }
    DLOG(INFO) << "Linearizing 2D image for writing in fits";
    linearize_2Dvector(image, &image1D);

    fits_movnam_hdu(fimg, IMAGE_HDU, (char*)extname.c_str(), 0, &status);
    errorMsg = "Error in moving to " + (string) extname + " extension in file ";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    
    fits_update_key(fimg, TINT, "NAXIS1", &naxes[0], NULL, &status);
    fits_update_key(fimg, TINT, "NAXIS2", &naxes[1], NULL, &status);
    
    DLOG(INFO) << "[DEBUG] Writing image with dimensions " << naxes[0] << " x " << naxes[1];
    fits_write_pix(fimg, DATA_TYPE, fpixel, imgSize, image1D.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing image data";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}


template <class T> int read_img(vector <vector <T> > &image, fitsfile *fimg, string extname, int DATA_TYPE){
    T* imageArray;
    vector <T> image1D;
    long ix=0, iy=0;
    int status=0;
    string errorMsg="";
    int maxdim=2;
    long naxes[2];
    long fpixel[2];
    long imgSize;
    fpixel[0] = fpixel[1]=1;

    fits_movnam_hdu(fimg, IMAGE_HDU, (char*)extname.c_str(), 0, &status);
    errorMsg = "Error in moving to " + (string) extname + " extension in file ";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    
    fits_get_img_size(fimg, maxdim, naxes, &status);
    if (status) {
        LOG(ERROR) << "Error in getting image size.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
    imgSize=naxes[0]*naxes[1];
    
    DLOG(INFO) << "[DEBUG] Reading image with dimensions " << naxes[0] << " x " << naxes[1];
    imageArray = new T [imgSize];
    fits_read_pix(fimg, DATA_TYPE, fpixel, imgSize, NULL, imageArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading image data";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    image1D.resize(naxes[0],0);
    image.resize(naxes[1], image1D);
    
    for(iy=0; iy<naxes[1]; iy++){
        for(ix=0; ix<naxes[0]; ix++){
            image[iy][ix] = imageArray[iy*naxes[0]+ix];
        }
    }
    
    return EXIT_SUCCESS;
}

//Write fits column
template <class T>
int write_fits_column(fitsfile *fptr, string colname, 
        int datatype, long firstrow, long firstelem, 
        vector <T> vec_data){
    int status=0;
    int colnum=0;
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fptr, datatype, colnum, firstrow, 1, vec_data.size(),
            vec_data.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing " << colname << "COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}

template <class T>
int write_fits_array_column(fitsfile *fptr, string colname, int datatype, long firstrow, long firstelem,
        vector < vector <T> > vec_data) {
    int status = 0;
    int colnum = 0;
    vector <T> vec_data1D;
    if (linearize_2Dvector(vec_data, &vec_data1D)) {
        LOG(ERROR) << "Error in linearizing 2D vector.";
        return EXIT_FAILURE;
    }
    
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_modify_vector_len(fptr, colnum, vec_data[0].size(), &status);
    if (status) {
        LOG(ERROR) << "Error in modifying vector length of " << colname << ".";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fptr, datatype, colnum, firstrow, firstelem, vec_data1D.size(), vec_data1D.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing " << colname << "column";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return status;
    
}
//Read fits column
/**
 * Reads and stores fits column values in a vector.
 * @param fptr
 * @param colname: Column to be read
 * @param datatype: Datatype of column
 * @param firstrow: start row number
 * @param firstelem: start element number
 * @param nelements: number of elements to be read. 
 *                   nelements=-1 => complete row will be read.
 * @param vec_data: vector to store column data.
 * @return 
 */
template <class T>
int read_fits_column(fitsfile *fptr, string colname, 
        int datatype, long firstrow, long firstelem, long nelements, 
        vector <T> &vec_data){
    int status=0;
    long i=0;
    long nrows=0;
    int colnum=0;
    T* temp_data;
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    if(nelements==-1){
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for " << colname << "COLUMN.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        nelements=nrows;
    }
    temp_data = new T [nelements];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading " << colname << "COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Assigning values to output vector
    for(i=0; i<nelements; i++){
        vec_data.push_back(temp_data[i]);
    }
    
    delete[] temp_data;
    return EXIT_SUCCESS;
}

template <class T>
int read_fits_columnN(fitsfile *fptr, string colname,
        int datatype, long firstrow, long firstelem, long nelements,
        vector <T> &vec_data) {
    int status = 0;
    long i = 0;
    long nrows = 0;
    int colnum = 0;
    ErrorHandler errHandler;
    T* temp_data;
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in getting column number for " + colname + "."; "";
        throw errHandler;
    }
    
    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            fits_read_errmsg(errHandler.fitsErrMsg);
            fits_get_errstatus(status, errHandler.fitsErrTxt);
            errHandler.fitsflag = true;
            errHandler.fitsErrorStatus = status;
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in getting number of rows for " + colname;
            throw errHandler;
        }
        nelements = nrows;
    }
    
    temp_data = new T [nelements];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading " + colname + "COLUMN.";
        throw errHandler;
    }

    //Assigning values to output vector
    for (i = 0; i < nelements; i++) {
        vec_data.push_back(temp_data[i]);
    }

    delete[] temp_data;
    return EXIT_SUCCESS;
}

//Read fits array column
template <class T>
int read_fits_array_column(fitsfile *fptr, string colname, int datatype, long firstrow, 
        long firstelem, long nelements,
        vector < vector <T> > &vec_data) {
    int status = 0;
    int colnum = 0;
    int typecode=0;
    long repeat =0;
    long width =0;
    long nrows = 0;
    long ivalue=0, ielement=0;
    vector <T> vec_data1D;
    T* tempData;


    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column type of " << colname << ".";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for " << colname << "COLUMN.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        nelements = nrows*repeat;
    }
    
    if(datatype!=typecode){
        LOG(ERROR) << "Improper datatype of column to be read.";
        return EXIT_FAILURE;
    }
    
    tempData = new T [nelements];
    
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            tempData, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading " << colname << "COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Assigning values to output vector
    for (ielement = 0; ielement < nelements; ielement+=repeat ) {
        for(ivalue=0; ivalue<repeat; ivalue++){
            vec_data1D.push_back(tempData[ielement + ivalue]);
        }
        vec_data.push_back(vec_data1D);
        vec_data1D.clear();
    }

    delete[] tempData;
    return EXIT_SUCCESS;
}

//Read fits array column new
template <class T>
int read_fits_array_columnN(fitsfile *fptr, string colname, int datatype, long firstrow, 
        long firstelem, long nelements,
        vector < vector <T> > &vec_data) {
    int status = 0;
    int colnum = 0;
    int typecode=0;
    long repeat =0;
    long width =0;
    ErrorHandler errHandler;
    long nrows = 0;
    long ivalue=0, ielement=0;
    vector <T> vec_data1D;
    T* tempData;

    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg ="Error in getting column number for " + colname;
        throw errHandler;
    }

    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in getting column type of " + colname + ".";
        throw errHandler;
    }


    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            fits_read_errmsg(errHandler.fitsErrMsg);
            fits_get_errstatus(status, errHandler.fitsErrTxt);
            errHandler.fitsflag = true;
            errHandler.fitsErrorStatus = status;
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in getting number of rows for " + colname + "COLUMN.";
            throw errHandler;
        }
        nelements = nrows*repeat;
    }
    
    if(datatype!=typecode){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT_FITS_DATATYPE;
        errHandler.errorMsg = "Improper datatype of column to be read.";
        throw errHandler;
    }
    
    tempData = new T [nelements];
    
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            tempData, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading " << colname << "COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Assigning values to output vector
    for (ielement = 0; ielement < nelements; ielement+=repeat ) {
        for(ivalue=0; ivalue<repeat; ivalue++){
            vec_data1D.push_back(tempData[ielement + ivalue]);
        }
        vec_data.push_back(vec_data1D);
        vec_data1D.clear();
    }

    delete[] tempData;
    return EXIT_SUCCESS;
}

template <class T>
int read_fits_column_if_available(fitsfile *fptr, string colname, int datatype,
         long firstrow, long firstelem, long nelements, vector <T> &vec_data) {
    int status = 0;
    long i = 0;
    long nrows = 0;
    int colnum = 0;
    string errorMsg="";
    char errMsg[100];
    T* temp_data;
    ErrorHandler errHandler;
    
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errWARNING;
        errHandler.fitsflag = true;
        errHandler.errorMsg = "Unable to read column " + colname + ".";
        throw errHandler;
    }
    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for " << colname << "COLUMN.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        nelements = nrows;
    }
    temp_data = new T [nelements];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading " << colname << "COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Assigning values to output vector
    for (i = 0; i < nelements; i++) {
        vec_data.push_back(temp_data[i]);
    }

    delete[] temp_data;
    return EXIT_SUCCESS; 
}

template <class T>
int read_fits_column_if_availableN(fitsfile *fptr, string colname, int datatype,
         long firstrow, long firstelem, long nelements, vector <T> &vec_data) {
    int status = 0;
    long i = 0;
    long nrows = 0;
    int colnum = 0;
    string errorMsg="";
    char errMsg[100];
    T* temp_data;
    ErrorHandler errHandler;
    
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errWARNING;
        errHandler.fitsflag = true;
        errHandler.errorMsg = "Unable to read column " + colname + ".";
        throw errHandler;
    }
    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            fits_read_errmsg(errHandler.fitsErrMsg);
            fits_get_errstatus(status, errHandler.fitsErrTxt);
            errHandler.fitsflag = true;
            errHandler.fitsErrorStatus = status;
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in getting number of rows for " + colname + "COLUMN.";
            throw errHandler;
        }
        nelements = nrows;
    }
    temp_data = new T [nelements];
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column " + colname;
        throw errHandler;
    }

    //Assigning values to output vector
    for (i = 0; i < nelements; i++) {
        vec_data.push_back(temp_data[i]);
    }

    delete[] temp_data;
    return EXIT_SUCCESS; 
}

int read_fits_string_column(fitsfile *fptr, string colname,
        int datatype, long firstrow, long firstelem, long nelements,
        vector <string> &vec_data);

int read_fits_string_columnN(fitsfile *fptr, string colname,
        int datatype, long firstrow, long firstelem, long nelements,
        vector <string> &vec_data);


//ERROR HANDLING

void logError(ErrorHandler errHandler);

//Statistics
bool max_pair_value(pair<float, long> pair1, pair<float, long> pair2);

template <class T>
void calculate_mode(vector<T> vecNumbers, T *mode, map<T,long> *modemap){
    ErrorHandler errHandler;
    vector<T> uniqueNumbers;
    typename std::vector<T>::iterator it;
    typename std::map<T,long>::iterator itmap;
    pair<T,long> maxPair;
    vector<T> modeValues; //all temperatures with maximum number of occurences.
    T sumModeValues=0.0;
    T meanModeValues=0.0;
    long iNumber=0;
    long nOccurences=0;
    
    //Error Handling
    if(vecNumbers.size()<=0) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Cannot calculate mode of n numbers when n=0.";
        throw errHandler;
    }
    //Getting number of unique values in the vector
    uniqueNumbers=vecNumbers;
    it=unique(uniqueNumbers.begin(), uniqueNumbers.end());
    uniqueNumbers.resize(distance(uniqueNumbers.begin(), it));
    
    for(iNumber=0; iNumber<uniqueNumbers.size(); iNumber++){
        nOccurences = count(vecNumbers.begin(), vecNumbers.end(), uniqueNumbers[iNumber]);
        (*modemap).insert(pair<T,long>(uniqueNumbers[iNumber], nOccurences));
    }
    
    
    //Getting maximum frequency pair #returns only a single pair, there is a
    //chance that many such pairs exist.
    maxPair = *max_element((*modemap).begin(), (*modemap).end(), max_pair_value);
    
    modeValues.clear();
    //Getting all the pairs with this maximum frequency
    for(itmap=(*modemap).begin(); itmap!=(*modemap).end(); ++itmap){
        if(itmap->second==maxPair.second){
            modeValues.push_back(itmap->first);
        }
    }
    
    sumModeValues = accumulate(modeValues.begin(), modeValues.end(), 0.0);
    meanModeValues = sumModeValues/modeValues.size();
    
    //returning mode Temperature
    *mode = meanModeValues;

}

template <class in_T, class out_T>
void assign_vector(vector < vector <in_T> > input, vector < vector <out_T> >  *output){
    long xsize=0;
    long ysize=0;
    long i=0, j=0;
    ErrorHandler errHandler;
    vector <out_T> tempOut;
    ysize = input.size();
    if(ysize==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Input vector size: 0";
        throw errHandler;
    } else {
    xsize = input[0].size();
    }
    
    tempOut.resize(xsize, 0.0);
    (*output).resize(ysize, tempOut);
    
    for(i=0; i<ysize; i++){
        for(j=0; j<xsize; j++){
            (*output)[i][j] = (out_T) input[i][j];
        }
    }
}

int get_quadsToProcessVec(string quadsToProcess, vector<int> *quadVec);

template <class T>
void print_vector(vector <T> vecData, int precision=10){
    int i=0;
    for(i=0; i<vecData.size(); i++){
        LOG(INFO) << setprecision(precision) << vecData[i];
    }
}
//FILE/DIRECTORY FUNCTIONS;
//*******************Function declarations*******************
/**
 * Function returns true if the File Exists
 * @param filename
 * @return 
 */
bool FileExists(char *filename);
//----------------------------------------------------------------

void files_exist(vector<string> filenames, vector<string> *nonExistentFiles);
/**
 * Function returns true if the directory exists
 * @param dirname
 * @return 
 */
bool DirExists(char *dirname);
//----------------------------------------------------------------
/**
 * Delete files if clobber=yes
 * @param filenames
 * @param clobber
 * @return 
 */
int deleteFile(vector<string> filenames, int clobber);
/**
 * Copies n keywords from fits file pointed by fptr1 to fptr2 
 * @param fptr1
 * @param fptr2
 * @param n
 * @param ...
 */
int getFiles(string dirName, vector<string> *filelist);
int getDirs(string dirName, vector<string> *subDirs);
int searchDir(string dirname, string pattern, vector<string> *fileList, vector<string> *dirList);

int getHeaderKey(const char* infile, int extno,const char*key,char* value);
vector <string> find_files_in_filelist(vector<string> filelist, string extname, vector <string> *extraChecks);

int updatelivetime(double gtstart, double gtstop,double *livetime,long tstarti, float fraction,bool add,long ntbins,double livetime_binsize);
int updatelivetime_gti(fitsfile *fgti,double *timearray, double *livetime,long ntbins,double livetime_binsize);
int re_updatelivetime(double gtstart, double gtstop,double *inlivetime,double *outlivetime,long tstarti, float fraction,bool add,long ntbins,double livetime_binsize);


#endif	/* COMMON_H */

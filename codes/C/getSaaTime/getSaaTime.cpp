//getSaaTime.cpp

// Function to compute SAA crossing times based on lat and lon values from
// CZTI MKF file.

// Mithun N P S (16/06/20) 

#include <stdio.h>
#include <fitsio.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

int getSaaTime(string mkfFile, vector <double> &saaCrossTimes);

int main()
{
    string mkfFile;
    int status,i;
    vector <double> saaCrossTimes;

    //mkfFile="/home/mithun/work/czti/CZTI_noise/data/9000000618/AS1G05_009T02_9000000618czt_level2.mkf";
    cout<<"Enter CZTI MKF file name: ";
    cin>>mkfFile;

    cout<<"Using MKF file: "<<mkfFile<<"\n";

    status=getSaaTime(mkfFile, saaCrossTimes);
    if(status)
    {
        cout<<"Error in getting SAA crossing times\n";
        return status;
    }

    cout<<"Number of SAA Crossings: "<<saaCrossTimes.size()<<"\n";
    for (i=0;i<saaCrossTimes.size();i++)
        printf("%lf\n",saaCrossTimes[i]);

    
    return EXIT_SUCCESS;    
}

int getSaaTime(string mkfFile, vector <double> &saaCrossTimes)
{

    int status=0,hdutype;
	int time_col,lat_col,lon_col;
	long nrows,i;
	double *time;
	float *lat,*lon;

    fitsfile *fout;

	// Read MKF file and get lat and lon with time

    fits_open_file(&fout, mkfFile.c_str(), READONLY, &status);
    if (status) {
        cout <<"***Error in opening file "<<mkfFile<<"\n";
        cout <<"***Error: Check the file format\n";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fout, 2, &hdutype, &status);
    if (status) {
         cout <<"***Error: Check the file format\n";
         cout <<"***Error in moving to data HDU in file "<<mkfFile<<"\n";
         fits_report_error(stderr, status);
         return (EXIT_FAILURE);
    }
    
    fits_get_colnum(fout,CASEINSEN,(char *)"TIME",&time_col,&status);
    fits_get_colnum(fout,CASEINSEN,(char *)"EARTHLAT",&lat_col,&status);
    fits_get_colnum(fout,CASEINSEN,(char *)"EARTHLON",&lon_col,&status);
    
    if (status) {
        cout <<"***Error: Check the file format\n";
        cout <<"***Error in getting Expected column numbers for file "<<mkfFile<<"\n";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_get_num_rows(fout, &nrows, &status);
    
    if(nrows==0){
        cout<<"***Error: No rows in input MKF file "<<mkfFile<<"\n";
        return EXIT_FAILURE;
    }
    
    time=(double *)malloc(sizeof(double)*nrows);
    lat=(float *)malloc(sizeof(float)*nrows);
    lon=(float *)malloc(sizeof(float)*nrows);
    
    fits_read_col(fout,TDOUBLE,time_col,1,1,nrows,NULL,time,NULL,&status);
    fits_read_col(fout,TFLOAT,lat_col,1,1,nrows,NULL,lat,NULL,&status);
    fits_read_col(fout,TFLOAT,lon_col,1,1,nrows,NULL,lon,NULL,&status);
    
    if(status)
    {
        cout<<"***Error in reading MKF file "<<mkfFile<<"\n";
        return EXIT_FAILURE;
    }
	
    fits_close_file(fout,&status);
    if (status) {
        cout <<"***Error in closing file "<<mkfFile<<"\n";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

	// Find SAA crossing times

    // SAA crossing is defined as the time when the spacecraft crosses 
    // the line lon=(((lat+6)*slope)-50.0) with slope=20.0/12.0 in the 
    // latitude, longitude space. Crossing is identified by evaluating 
    // the line equation with two successive lat,lon and when it changes
    // sign from negative to positive.

    saaCrossTimes.clear();

    float slope=20.0/12.0;
    float lval1=0,lval2=0;

    lval1=(lon[0]-(((lat[0]+6)*slope)-50.0));

    for(i=1;i<nrows-1;i++)
    {
        lval2=(lon[i]-(((lat[i]+6)*slope)-50.0));
        if(lval1<0 && lval2 >=0)
            saaCrossTimes.push_back(time[i]);
        lval1=lval2;    
    }

    free(time);
    free(lat);
    free(lon);

    return EXIT_SUCCESS;    
}

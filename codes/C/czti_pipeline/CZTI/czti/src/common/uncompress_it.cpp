#include <iostream>
#include <fstream>
#include <sstream>
#include <fitsio.h>
#include<string.h>

#define NUMQUAD 4
#define NUMROWS 8100
#define NUMCOLS 8100
#define COMPCOLS 270
int restore_mask_char(char *compmask_filename, unsigned char **mask);
int uncompress(int *maskpattern,unsigned char *uncompressed_field, int i);
void writeFits(char *filename,unsigned char** pixels);

using namespace std;


int main(int argc, char *argv[]) 
{
	int i;	
	unsigned char **mask;

  mask=new unsigned char*[NUMQUAD];
  for(int i=0;i<NUMQUAD;i++){
  	mask[i]=new unsigned char[NUMROWS*NUMCOLS];
	}
	restore_mask_char(argv[1],mask);
	return 0;
}

int restore_mask_char(char *compmask_filename, unsigned char **mask)
{
	int status=0, hdutype, *buffer_pix=NULL, anynull, i,j,k;
	long naxes[2], fpixel=1, nbuffer, npixels=NUMROWS*COMPCOLS, f_pix[2]={1,1};
	float nullval=0;
	unsigned char *uncompressed_field; // TO STORE 30 BITS EACH AS A SINGLE CHARACTER
	//unsigned char *mask[4];

  for(int i=0;i<NUMQUAD;i++){
   for(j=0;j<(NUMROWS*NUMCOLS);j++){
    mask[i][j]='0';
   }
  }
	fitsfile *fptr;
	buffer_pix = (int*)malloc(sizeof(int)*NUMROWS*COMPCOLS);
	uncompressed_field = (unsigned char*)malloc(sizeof(char)*30);
	//for(i=0;i<4;i++) {mask[i]=(unsigned char*)malloc(sizeof(char)*NUMROWS*NUMCOLS);}
			
	fits_open_file(&fptr, compmask_filename, READONLY, &status);
	for(k=2;k<6;k++){
		fits_movabs_hdu(fptr, k, &hdutype, &status);
		fits_read_pix(fptr, TINT, f_pix, npixels, &nullval, buffer_pix, &anynull, &status);
		for(i=0;i<npixels;i++){
			uncompress(buffer_pix, uncompressed_field, i);
			for(j=0; j<30; j++){
				mask[k-2][i*30+j]=uncompressed_field[j];		
			}		
		}
	}
	fits_close_file(fptr, &status);
	
	writeFits("!check.fits",mask);
	return 0;
}

int uncompress(int *maskpattern,unsigned char *uncompressed_field, int i)
{
	int k=0,val=0, temp_val=0;
	//cout << maskpattern[i] << endl;
	for(k=0; k<30; k++){
		temp_val = 1<<k;
		val= maskpattern[i]&temp_val;
		if(val){
		uncompressed_field[k]= '1';
		}
		else uncompressed_field[k]='0';
	}	
	return 0;	
}

void writeFits(char *filename,unsigned char** pixels)
{
	fitsfile *fptr;
	int bitpix   =  BYTE_IMG;
    	long naxis    =   2;  
    	int status=0,i=0;         
	long fpixel[2] = {1,1};                      
	long naxes[2] = { 8100,8100 };
	long nelements = naxes[0] * naxes[1];  
	long zero=0;
	char command[100],temp[100],extensionName[100];

	fits_create_file(&fptr, filename, &status);
	//fits_create_img(fptr,bitpix, naxis, naxes, &status);
	if(status){
		cerr<<"\nCouldn't create file "<<filename; exit(1);
	}
	for(int i=0;i<4;i++){
		fits_create_img(fptr,bitpix, naxis, naxes, &status);
		if(status){
			cerr<<"\nCouldn't create image at loop number "<<i;  exit(1);
		}
		sprintf(temp,"%d",i);
		strcpy(extensionName,"Q");
		strcat(extensionName,temp);
		fits_write_key(fptr, TSTRING, "EXTNAME", extensionName, "", &status); 
		if(status){
			cerr<<"\nCouldn't write EXTNAME at loop number "<<i; exit(1);
		}
		fits_write_pix(fptr, TBYTE, fpixel, nelements, pixels[i], &status);
	}
	fits_close_file(fptr,&status);
}


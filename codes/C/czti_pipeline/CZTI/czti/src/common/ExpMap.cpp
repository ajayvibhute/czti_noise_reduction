
#include"ExpMap.h"
#include<iomanip>

using namespace std;
/*****************************************************************************
//Is is assumed that fits file for Mask will have four extensions for each quadrant with 
 * 8175 x 325 pixels containing mask pattern in compressed form
 *****************************************************************************/
//function definitions for member functions of 'Mask' class
 
Mask::Mask(){
    strcpy(maskfile,"\0");
}

int Mask::read(char *maskfile){
    LOG(INFO)<<"Inside Mask::read()";
    strcpy(this->maskfile,maskfile); 
    string maskFile;
    mask=new unsigned char *[NUMQUAD];
    for(int i=0;i<NUMQUAD;i++){
        mask[i]=new unsigned char[QUAD_MASKSIZE*QUAD_MASKSIZE];
        if(mask[i]==NULL){
            LOG(ERROR)<<"***Out of Memory Error***";
            LOG(ERROR)<<"***Memory allocation failed for storing mask file***";
            return (EXIT_FAILURE);
        }
    }
    
    //---------------------------------------------------------------
     int status=0, hdutype, *buffer_pix=NULL, anynull, i,j,k;
     long naxes[2], fpixel=1, nbuffer, npixels=COMPROWS*COMPCOLS, f_pix[2]={1,1};
     float nullval=0;
     unsigned char *uncompressed_field; // TO STORE 30 BITS EACH AS A SINGLE CHARACTER
        //unsigned char *mask[4];

    for(i=0;i<NUMQUAD;i++){
        for(j=0; j<(QUAD_MASKSIZE*QUAD_MASKSIZE);j++)
                mask[i][j]=0;
    }
   
    fitsfile *fptr;
    buffer_pix = (int*)malloc(sizeof(int)*COMPROWS*COMPCOLS);
    uncompressed_field = (unsigned char*)malloc(sizeof(char)*25);
    //for(i=0;i<4;i++) {mask[i]=(unsigned char*)malloc(sizeof(char)*NUMROWS*NUMCOLS);}
    LOG(INFO) << "STATUS IS :" << status;

    maskFile = (string) maskfile;

    fits_open_file(&fptr, (char*) maskFile.c_str(), READONLY, &status);
    if(status){
        LOG(ERROR)<<"***Error opening file "<<maskFile<<" ***";
        return (EXIT_FAILURE);
    }
    for(k=2;k<6;k++){
            LOG(INFO)<<"Reading Quadrant "<<k-1<<"............";   
            fits_movabs_hdu(fptr, k, &hdutype, &status);
            if(status){
                LOG(ERROR)<<"***Error moving to HDU  "<<k<<"  in file "<<maskFile<<" ***";
                return (EXIT_FAILURE);
            }
            fits_read_pix(fptr, TINT, f_pix, npixels, &nullval, buffer_pix, &anynull, &status);
            if(status){
                  LOG(ERROR)<<"***Error reading pixels in file "<<maskFile<<" ***";
                  return (EXIT_FAILURE);
            }
            for(i=0;i<npixels;i++){
                    uncompress(buffer_pix, uncompressed_field, i);
                    for(j=0; j<COMP_FIELD_SIZE; j++){
                            mask[k-2][i*COMP_FIELD_SIZE+j]=uncompressed_field[j];
                            //cout<<"\n"<<mask[k-1][i*30+j];
                    }
            }
    }
    /*-------------------------------------------------------------*/
    /********Print uncompressed mask in mask file */
    /*-------------------------------------------------------------*/
    fits_close_file(fptr, &status);
//    char filename[12];
//     sprintf(filename,"mask_check.fits");
//    write(filename);
    /*-----------------------------------------------------*/
    
    
    if(status){
                LOG(ERROR)<<"***Error closing file "<<maskFile<<" ***";
                return (EXIT_FAILURE);
        }
    free(buffer_pix);
    free(uncompressed_field);
     //code for reading mask file here
    return (EXIT_SUCCESS);
}

int Mask::uncompress(int *maskpattern,unsigned char *uncompressed_field, int i)
{
        int k=0,val=0, temp_val=0;
        //cout << maskpattern[i] << endl;
        for(k=0; k<COMP_FIELD_SIZE; k++){
                temp_val = 1<<k;
                val= maskpattern[i]&temp_val;
                if(val){
                uncompressed_field[k]= 1;
                }
                else uncompressed_field[k]=0;
        }
        return (EXIT_SUCCESS);
}

int Mask::write(char *filename)
{
        fitsfile *fptr;
        int bitpix   =  BYTE_IMG;
        long naxis    =   2;
        int status=0,i=0;
        long fpixel[2] = {1,1};
        long naxes[2] = { 8175,8175 };
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
                fits_write_pix(fptr, TBYTE, fpixel, nelements, mask[i], &status);
        }
        fits_close_file(fptr,&status);
        return (EXIT_SUCCESS);
}


/* getDetMask below is taking quadrant id , detector id as the input and the producing a 2D array of detmask.*/
int Mask::getDetMask(int quadid, int detid,float **detmask){
    
    //assuming mask values are stored in mask array
    int x1,y1,i,j,x,y;
    int row,col;           //row,col of detector in the quadrant
    if(detid<0){
        LOG(ERROR)<<"***"<<detid<<" :Detector id cannot be negative***";
        return (EXIT_FAILURE);
    }
    col=detid%NUM_DET_PER_QUAD_PER_ROW;
    row=detid/NUM_DET_PER_QUAD_PER_ROW;

    y1=row*(DET_MASKSIZE+NPIXELS_GAP);
    x1=col*(DET_MASKSIZE+NPIXELS_GAP);
    
    //cout<<"x1:"<<x1<<"  y1:"<<y1<<"  posx:"<<posx<<"  posy:"<<posy;

    int index=0;
    for(i=y1,y=0;y<DET_MASKSIZE;i++,y++){
        for(j=x1,x=0;x<DET_MASKSIZE;j++,x++){
            index=i*QUAD_MASKSIZE+j;
            detmask[y][x]=(float)mask[quadid][index];
        }
    }
    return (EXIT_SUCCESS); 
}

Mask::~Mask(){
    for(int i=0;i<NUMQUAD;i++)
        delete[] mask[i];
}

//function definitions for MaskWeighting functions

 MaskWeighting::MaskWeighting(char *maskfile,double thetax,double thetay,float flux){
     int i,  j;
     openf=new float[XSIZE*YSIZE];
      
     //initializing variable to store openfraction quadrant wise
    
    openf_quad=new float *[NUMQUAD];
    if(openf_quad==NULL)  { LOG(ERROR)<<"Out of memory error"; exit (EXIT_FAILURE); }
    for(int i=0;i<NUMQUAD;i++){
        openf_quad[i]=new float[NUMQUAD,PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD];
        if(openf_quad[i]==NULL)   { LOG(ERROR)<<"Out of memory error"; exit (EXIT_FAILURE); }
    } 

      
     for (i=0;i<NUMQUAD; i++){
         for (j=0; j< PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD; j++)
         openf_quad[i][j]=0.0;
     }
     // open fraction quadrant wise array initialized.

     
     area=new float[XSIZE*YSIZE];
     if(area==NULL || openf==NULL){
         LOG(ERROR)<<"***Out of memory Error***";
         exit(EXIT_FAILURE);
     }
     strcpy(this->maskfile,maskfile);
     sourceflux=flux;
     this->thetax=thetax;
     this->thetay=thetay;
     //computing deltax and deltay
    dx=MASKHEIGHT*tan(thetax);
    dy=MASKHEIGHT*tan(thetay);
    LOG(INFO)<<"dx:"<<dx<<"   dy:"<<dy;
   
 }
 
 int MaskWeighting::run(){
     LOG(INFO)<<"Computing Area ...........";
     computeArea();
     
     LOG(INFO)<<"Reading Maskfile ...........";
     if(read(maskfile)){
          LOG(ERROR)<<"***Error reading maskfile***";
         return (EXIT_FAILURE);
     }
     
     LOG(INFO)<<"Computing Open fraction using mask weighting ...........";
     if(computeOpenFraction()){
          LOG(ERROR)<<"***Error Computing Open fraction***";
         return (EXIT_FAILURE);
     }
     return (EXIT_SUCCESS);
}
 
 MaskWeighting::~MaskWeighting(){
     delete[]  openf, area;
 }
 
 
 
 
// computeOpenFraction()- function to compute the open fraction of the all quadrants and store it in an array openfraction[4][64x64]
 
 int MaskWeighting::computeOpenFraction(){
     
     // initializing detector mask and detector_open_fraction arrays to store uncompressed mask for each detector and its openfraction.
     
     float **detmask=allocateMemory<float>(DET_MASKSIZE,DET_MASKSIZE); //stores detector mask as obtained be getDetMask
     float **detopenfraction=allocateMemory<float>(PIXELS_PER_ROW_DET,PIXELS_PER_COL_DET);
     LOG(INFO) << "detmask and detopenfraction arrays created.";
     //float **openf_quad=allocateMemory<float>(NUMQUAD, PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD);
     if(detmask==NULL || detopenfraction==NULL){
         LOG(ERROR)<<"***Out of memory Error***";
         return (EXIT_FAILURE);
     }
     //Initialization of detmask and detopenfraction completed
     /****************************************************************************************/
     
     int qid,detid,pixid;
     int detrowno,detcolno;
     int row,col,index;
      char temp[512];
     for(qid=0;qid<NUMQUAD;qid++){
         for(int detid=0;detid<NUM_DET_PER_QUAD;detid++){
              getDetMask(qid,detid,detmask);
              shiftMaskDet(detmask);         //shifting the mask of detector with dx,dy
              getOpenFractionDet(detmask,detopenfraction);        //computing open fraction with shifted mask

              detrowno=detid/NUM_DET_PER_QUAD_PER_ROW;
              detcolno=detid%NUM_DET_PER_QUAD_PER_ROW;
              row=detrowno*PIXELS_PER_ROW_DET;
              col=detcolno*PIXELS_PER_COL_DET;
              
              for(int pixrowindex=0,x=row;pixrowindex<PIXELS_PER_COL_DET;pixrowindex++,x++){
                    for(int pixcolindex=0,y=col;pixcolindex<PIXELS_PER_ROW_DET;pixcolindex++,y++){

                        index=x*PIXELS_PER_ROW_QUAD+y;
                        openf_quad[qid][index]=detopenfraction[pixrowindex][pixcolindex];
                                
                    }
               }
         }

     }
      makeFullImage<float>(openf_quad[0],openf_quad[1],
                 openf_quad[2],openf_quad[3],openf);
      freeMemory(detmask,DET_MASKSIZE);
      freeMemory(detopenfraction,PIXELS_PER_ROW_DET);
      return (EXIT_SUCCESS);
 }
 
  int MaskWeighting::getMaskWeight(double *maskweight,int size){
      
     for(int i=0;i<size;i++){
         if(openf[i]==NAN){
             LOG(ERROR)<<"***Open Fraction is NAN at ("<<size/XSIZE<<", "<<size%YSIZE<<")***";
             return (EXIT_SUCCESS);
         }
         maskweight[i]=2*(double)openf[i];
     }
     return (EXIT_SUCCESS);
 }
  
  int MaskWeighting::getMaskWeightQuadWise(double **mskwt_quad, int size){

      for(int i=0; i<NUMQUAD; i++){
          for(int j=0; j<size; j++){
             if(openf_quad[i][j]==NAN){
                LOG(ERROR)<<"***Open Fraction is NAN at ("<<size/XSIZE<<", "<<size%YSIZE<<")***";
                return (EXIT_SUCCESS);
                }
              mskwt_quad[i][j]=2*(double)openf_quad[i][j] ;
          }
//          char filename[256];
//          sprintf(filename,"!file%d.fits",i);
//          writeImgD(filename,  mskwt_quad[i],64,64);  
      }
      
      return(EXIT_SUCCESS);
  }
 
 
 int MaskWeighting::getShadow(float *shadow,int size){
     
    if(dx>DETECTORWIDTH || dy>DETECTORWIDTH){
        LOG(ERROR)<<"***Shift greater than the limit***";
        return (EXIT_FAILURE);
    }
       
    double C=sourceflux/totalarea;
    LOG(INFO)<<"Constant C value:"<<C;
    double sourceangle=atan(sqrt(tan(thetax)*tan(thetax)+tan(thetay)*tan(thetay)));
    LOG(INFO)<<"Source Angle:"<<sourceangle;
    //cout<<"shadow";
    for(int i=0;i<YSIZE;i++){
        for(int j=0;j<XSIZE;j++){
            shadow[i*YSIZE+j]=C*cos(sourceangle)*openf[i*YSIZE+j]*area[i*YSIZE+j];
            if(shadow[i*YSIZE+j]==NAN){
                LOG(ERROR)<<"***Error in generating shadow***";
                return (EXIT_FAILURE);
            }
            //cout<<"|"<<shadow[i][j];
        }
    }
    return (EXIT_SUCCESS);
 }
 
 int MaskWeighting::getExposureMap(float* expmap, int size){
          
     memcpy(expmap,openf,XSIZE*YSIZE*sizeof(float));          //Open Fraction and Exposure Map are same
     return (EXIT_SUCCESS);
 }
 
   int MaskWeighting::getExposureMapQuadWise(float **openfraction_quad){

      for(int i=0; i<NUMQUAD; i++){
          for(int j=0; j<PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD; j++){
              openfraction_quad[i][j]=openf_quad[i][j];
          }
      }
      return(EXIT_SUCCESS);
  }
 
 
//function to shift the detector mask array by dx and dy
//int shiftMaskDet(float mask[DET_MASKSIZE][DET_MASKSIZE],double delta_x,double delta_y){
int MaskWeighting::shiftMaskDet(float **mask){  
//cout<<"Inside shiftMaskDet() - ExpMap.cpp";   
    double delta_x,delta_y;
    
    if(dx==0 && dy==0){
        return (EXIT_SUCCESS);
    }
    
    delta_x=dx,delta_y=dy;                 //converting to micrometer   
    int shiftx_pix=(delta_x)<0 ? (delta_x*(-1)/PIX_SIZ_MASK) : (delta_x/PIX_SIZ_MASK);  //shift in pixels
    int shifty_pix=(delta_y)<0 ? (delta_y*(-1)/PIX_SIZ_MASK) : (delta_y/PIX_SIZ_MASK);  //shift in pixels
                               //for one detector
    
    //LOG(INFO)<<"Shift in pixels : "<<shiftx_pix<<" , "<<shifty_pix;
    int i,j,x,y;
    float **shiftedmask=new float*[DET_MASKSIZE];
    if(shiftedmask==NULL){
        LOG(ERROR)<<"***Out of Memory**";
        return (EXIT_FAILURE);
    }
    for(i=0;i<DET_MASKSIZE;i++){
        shiftedmask[i]=new float[DET_MASKSIZE];
        if(shiftedmask[i]==NULL){
        LOG(ERROR)<<"***Out of Memory***";
        return (EXIT_FAILURE);
     }
    }
    //float shiftedmask[DET_MASKSIZE][DET_MASKSIZE];
    for(i=0;i<DET_MASKSIZE;i++){
        for(j=0;j<DET_MASKSIZE;j++){
           shiftedmask[i][j]=0;
        }
    }

    int rowstart,rowend,colstart,colend;  //for shifted mask
    int xstart,ystart,xend,yend;
    if(delta_x>0){
        colstart=0;
        colend=DET_MASKSIZE-shiftx_pix;
        xstart=shiftx_pix;
    }
    else{
        colstart=shiftx_pix;;
        colend=DET_MASKSIZE;;
        xstart=0;
    }

    if(delta_y>0){
        rowstart=0;
        rowend=DET_MASKSIZE-shifty_pix;
        ystart=shifty_pix;
    }
    else{
        rowstart=shifty_pix;
        rowend=DET_MASKSIZE;
        ystart=0;
       
    }

    
  for(i=rowstart,y=ystart;i<rowend;i++,y++){
        for(j=colstart,x=xstart;j<colend;j++,x++){
            shiftedmask[i][j]=mask[y][x];
            if(i>DET_MASKSIZE || j>DET_MASKSIZE){
                LOG(ERROR)<<"Index out of range for shifted mask";
               return(EXIT_FAILURE);
            }
             if(y>DET_MASKSIZE || x>DET_MASKSIZE){
                LOG(ERROR)<<"Index out of range for  mask";
                return(EXIT_FAILURE);
            }
            
        }
    }
    //copying the shifted mask back to original array
    //cout<<"Shifted mask";
    for(i=0;i<DET_MASKSIZE;i++){
        for(j=0;j<DET_MASKSIZE;j++){
            mask[i][j]=shiftedmask[i][j];
            //cout<<mask[i][j]<<endl;
        }
    }
    
    //LOG(INFO)<<"Copied shifted mask to original mask";
    
    for(i=0;i<DET_MASKSIZE;i++)  delete[] shiftedmask[i];
    delete[] shiftedmask;
    
    return (EXIT_SUCCESS);
}
    
//int getOpenFraction(float shiftedmask[DET_MASKSIZE][DET_MASKSIZE],
//        float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET]){


int MaskWeighting::getOpenFractionDet(float **shiftedmask,float **openFraction){
     //cout<<"Inside getOpenFraction() - ExpMap.cpp";   
     int i,j,closecount=0,x,y;
     double totpix=0;
     float openf=0;
     int xstart,xend,ystart,yend;
     for(i=0;i<PIXELS_PER_COL_DET;i++){
        for(j=0;j<PIXELS_PER_ROW_DET;j++){
           closecount=0;
           if(i==0) { xstart=0;  xend=114; }
           else if(i==15)  { xstart=114+14*123;  xend=228+14*123; }
           else { xstart=114+(i-1)*123;   xend=114+(i-1)*123+123; }

           if(j==0)  { ystart=0;   yend=114; }
           else if(j==15) { ystart=114+14*123;  yend=228+14*123; }
           else { ystart=114+(j-1)*123;  yend=114+(j-1)*123+123; }

           if((i==0 && j==0)||(i==15 && j==0)||(i==0 && j==15)||(i==15 && j==15))  totpix=114.0*114.0;
           else if((i==0 && 0<j<15)||(i==15 && 0<j<15)||(j==0 && 0<i<15)||(j==15 && 0<i<15))  totpix=114.0*123.0;
           else totpix=123.0*123.0;

           for(x=xstart;x<xend;x++){
                for(y=ystart;y<yend;y++){
                       if(shiftedmask[x][y]==CLOSE)
                                closecount++;
                }
           }
        openFraction[i][j]=1.0-((float)closecount/(float)totpix);
        }
    }
    return (EXIT_SUCCESS);
}


void MaskWeighting::computeArea(){
    //cout<<"Inside getArea() - ExpMap.cpp";               
    int i,j;
    float l,b; 
    for(i=0;i<YSIZE;i++){
        for(j=0;j<XSIZE;j++){
            switch(i){
                case 0:   
                case 15:   
                case 16: 
                case 31:
                case 32: 
                case 47:
                case 48:
                case 63: 
                case 64:
                case 79:
                case 80:
                case 95:
                case 96:
                case 111:
                case 112:
                case 127:   l=SMALL_PIXEL_WIDTH; //in mm
                break;
                default: l=BIG_PIXEL_WIDTH;  
            }
            switch(j){
                case 0:
                case 15:
                case 16:
                case 31:
                case 32:
                case 47:
                case 48:
                case 63:
                case 64:
                case 79:
                case 80:
                case 95:
                case 96:
                case 111:
                case 112:
                case 127:   b=SMALL_PIXEL_WIDTH;     //in mm
                break;
                default:  b=BIG_PIXEL_WIDTH;
            }
            area[i*YSIZE+j]=l*b;
        }
    } 
    totalarea=DETECTORWIDTH*DETECTORWIDTH*TOTAL_DET_MODULES;                    //in sq mm
}

int allToQuad(float *value, float **value_quad){
    int index0=0, index1=0, index2=0, index3=0;
    for(int i=0; i<64; i++){
        for(int j=0; j<64;j++){
            value_quad[0][index0]=value[i*128+j];
            index0++;
        }
    }   
    
    for(int i=0; i<64; i++){
        for(int j=64; j<128; j++){
            value_quad[1][index1]= value[i*128+j];
            index1++;
        }
    }
    
    for(int i=64; i<128; i++){
        for(int j=64; j<128; j++){
            value_quad[2][index2]= value[i*128+j];
            index2++;
        }
    }
    
    for(int i=64; i<128; i++){
        for(int j=0; j<64; j++){
            value_quad[3][index3]= value[i*128+j];
            index3++;
        }
    }    
}

int findKnownSources(char *catalog_filename, Mkf mkf, Teldef teldef, int extnumCatalog, 
        char* RA_colname, char* DEC_colname, char* FLUX_colname,
        vector<float> &RA, vector<float> &DEC, vector <float> &FLUX, 
        vector <float> &thetax, vector <float> &thetay) {

    int i = 0; //counter variable
    int status = 0; //status variable
    long nrows = 0;
    int hdutype = 0; // to store type of HDU
    float *flux, *ra, *dec, *theta_x, *theta_y; //to store RA, DEC, flux values from catalog file.
    int ra_colnum, dec_colnum, flux_colnum; // to store RA, DEC and flux column number from catalog file.

    /* temporary variables to store:
            1. ra_radian: to store RA value in radians
            2. dec_radian: to store DEC value in radians
            3. thetax_radian: to store camera coordinate theta_x in radian
            4. thetay_radian: to store camera coordinate theta_y in radian
            5. dx: shift in x direction
            6. dy: shift in y direction
     */
    float ra_radian, dec_radian, thetax_radian, thetay_radian, dx, dy;
    float nX, nY, nZ = 0.0; //to store inertial vector
    float detX, detY, detZ=0.0; //to store CZTI body vector
    
    fitsfile *fcatalog; //fits file pointer to catalog file

    fits_open_file(&fcatalog, catalog_filename, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error opening catalog file ***";
        return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fcatalog, extnumCatalog, &hdutype, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error while moving to required HDU number ***";
        return (EXIT_FAILURE);
    }
    if (hdutype != BINARY_TBL) {
        LOG(ERROR) << "*** Extension number " << extnumCatalog << " of " << catalog_filename << "is not a binary table. ***";
        return (EXIT_FAILURE);
    }

    fits_get_num_rows(fcatalog, &nrows, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error while getting number of rows ***";
        return (EXIT_FAILURE);
    }
    LOG(INFO) << "Total Sources in Catalog:" << nrows;

    flux = new float[nrows];
    ra = new float[nrows];
    dec = new float[nrows];
    theta_x = new float[nrows];
    theta_y = new float[nrows];

    // READING DATA FROM CATALOG
    fits_get_colnum(fcatalog, CASEINSEN, RA_colname, &ra_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fcatalog, CASEINSEN, DEC_colname, &dec_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fcatalog, CASEINSEN, FLUX_colname, &flux_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_read_col(fcatalog, TFLOAT, ra_colnum, 1, 1, nrows, NULL, ra, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcatalog, TFLOAT, dec_colnum, 1, 1, nrows, NULL, dec, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcatalog, TFLOAT, flux_colnum, 1, 1, nrows, NULL, flux, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //CATALOG DATA READ.

    fits_close_file(fcatalog, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error closing " << catalog_filename << " ***";
        return (EXIT_FAILURE);
    }

    // finding sources in the field of view
    for (i = 0; i < nrows; i++) {
        ra_radian = ra[i]*M_PI/180;
        dec_radian = dec[i]*M_PI/180;
        
        //getting inertial vectors
        if(to_nX_nY_nZ(ra_radian, dec_radian, nX, nY, nZ)){
            LOG(ERROR) << "Error in converting from inertial vectors nX, nY & nZ to Ra & Dec values.";
            return EXIT_FAILURE;
        }
        
        //transforming inertial vectors into body vectors
        if(get_body_vector(mkf, teldef, nX, nY, nZ, detX, detY, detZ)){
            LOG(ERROR) << "Error in transforming inertial vector into body vector.";
            return EXIT_FAILURE;
        }
        
        //getting thetax_radian & thetay_radian
        if(to_thetaX_thetaY(detX, detY, detZ, thetax_radian, thetay_radian)){
            LOG(ERROR) << "Error in converting from CZTI body vector to corresponding thetaX and thetaY.";
            return EXIT_FAILURE;
        }
        
        
        theta_x[i] = thetax_radian;
        theta_y[i] = thetay_radian;
        dx = (-1) * MASKHEIGHT * tan(thetax_radian);
        dy = (-1) * MASKHEIGHT * tan(thetay_radian);
        //cout<<"Shift of Mask  --- dx:"<<dx<<"  dy:"<<dy;
        if (abs(dx) < DETECTORWIDTH && abs(dy) < DETECTORWIDTH) {
            LOG(INFO) << "Catalog number of source found to be in field of view:" << i+1; //debug
            RA.push_back(ra_radian);
            DEC.push_back(dec_radian);
            FLUX.push_back(flux[i]);
            thetax.push_back(theta_x[i]);
            thetay.push_back(theta_y[i]);
        }

    }

    return EXIT_SUCCESS;

}

int findKnownSources(char *catalog_filename, string aspectFileName, int extnumCatalog,
        char* RA_colname, char* DEC_colname, char* FLUX_colname,
        vector<float> &RA, vector<float> &DEC, vector <float> &FLUX,
        vector <float> &thetax, vector <float> &thetay) {

    int i = 0; //counter variable
    int status = 0; //status variable
    long nrows = 0;
    int hdutype = 0; // to store type of HDU
    float *flux, *ra, *dec, *theta_x, *theta_y; //to store RA, DEC, flux values from catalog file.
    int ra_colnum, dec_colnum, flux_colnum; // to store RA, DEC and flux column number from catalog file.

    /* temporary variables to store:
            1. ra_radian: to store RA value in radians
            2. dec_radian: to store DEC value in radians
            3. thetax_radian: to store camera coordinate theta_x in radian
            4. thetay_radian: to store camera coordinate theta_y in radian
            5. dx: shift in x direction
            6. dy: shift in y direction
     */
    float ra_radian, dec_radian, thetax_radian, thetay_radian, dx, dy;
    float nX, nY, nZ = 0.0; //to store inertial vector
    float detX, detY, detZ = 0.0; //to store CZTI body vector

    fitsfile *fcatalog; //fits file pointer to catalog file

    fits_open_file(&fcatalog, catalog_filename, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error opening catalog file ***";
        return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fcatalog, extnumCatalog, &hdutype, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error while moving to required HDU number ***";
        return (EXIT_FAILURE);
    }
    if (hdutype != BINARY_TBL) {
        LOG(ERROR) << "*** Extension number " << extnumCatalog << " of " << catalog_filename << "is not a binary table. ***";
        return (EXIT_FAILURE);
    }

    fits_get_num_rows(fcatalog, &nrows, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error while getting number of rows ***";
        return (EXIT_FAILURE);
    }
    LOG(INFO) << "Total Sources in Catalog:" << nrows;

    flux = new float[nrows];
    ra = new float[nrows];
    dec = new float[nrows];
    theta_x = new float[nrows];
    theta_y = new float[nrows];

    // READING DATA FROM CATALOG
    fits_get_colnum(fcatalog, CASEINSEN, RA_colname, &ra_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fcatalog, CASEINSEN, DEC_colname, &dec_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fcatalog, CASEINSEN, FLUX_colname, &flux_colnum, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_read_col(fcatalog, TFLOAT, ra_colnum, 1, 1, nrows, NULL, ra, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcatalog, TFLOAT, dec_colnum, 1, 1, nrows, NULL, dec, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcatalog, TFLOAT, flux_colnum, 1, 1, nrows, NULL, flux, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //CATALOG DATA READ.

    fits_close_file(fcatalog, &status);
    if (status) {
        fits_report_error(stderr, status);
        LOG(ERROR) << "*** Error closing " << catalog_filename << " ***";
        return (EXIT_FAILURE);
    }

    // finding sources in the field of view
    for (i = 0; i < nrows; i++) {
        ra_radian = ra[i] * M_PI / 180;
        dec_radian = dec[i] * M_PI / 180;
//        //getting inertial vectors
//        if (to_nX_nY_nZ(ra_radian, dec_radian, nX, nY, nZ)) {
//            LOG(ERROR) << "Error in converting from inertial vectors nX, nY & nZ to Ra & Dec values.";
//            return EXIT_FAILURE;
//        }
//
//        //transforming inertial vectors into body vectors
//        if (get_body_vector(mkf, teldef, nX, nY, nZ, detX, detY, detZ)) {
//            LOG(ERROR) << "Error in transforming inertial vector into body vector.";
//            return EXIT_FAILURE;
//        }
//
//        //getting thetax_radian & thetay_radian
        if (to_thetaX_thetaY(aspectFileName, ra_radian, dec_radian, thetax_radian, thetay_radian)) {
            LOG(ERROR) << "Error in converting from CZTI body vector to corresponding thetaX and thetaY.";
            return EXIT_FAILURE;
        }


        theta_x[i] = thetax_radian;
        theta_y[i] = thetay_radian;
        dx = (-1) * MASKHEIGHT * tan(thetax_radian);
        dy = (-1) * MASKHEIGHT * tan(thetay_radian);
        //cout<<"Shift of Mask  --- dx:"<<dx<<"  dy:"<<dy;
        if (abs(dx) < DETECTORWIDTH && abs(dy) < DETECTORWIDTH) {
            LOG(INFO) << "Catalog number of source found to be in field of view:" << i + 1; //debug
            RA.push_back(ra_radian);
            DEC.push_back(dec_radian);
            FLUX.push_back(flux[i]);
            thetax.push_back(theta_x[i]);
            thetay.push_back(theta_y[i]);
        }

    }

    return EXIT_SUCCESS;

}



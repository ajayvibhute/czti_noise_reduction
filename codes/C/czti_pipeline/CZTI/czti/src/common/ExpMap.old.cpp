
#include"ExpMap.h"
#include<iomanip>


using namespace std;
/*****************************************************************************
//Is is assumed that fits file for Mask will have 4 extensions for 4 quadrants
 * The masks for 4 quadrants would be in hdu 2, 3, 4 and 5
 *****************************************************************************/
//function definitions for member functions of 'Mask' class
 
Mask::Mask(){
    strcpy(maskfile,"\0");
    maskquad=new float*[QUAD_MASKSIZE];
    if(maskquad==NULL){
        LOG(ERROR)<<"***Out of Memory Error***";
        exit(1);
    }
    for(int i=0;i<QUAD_MASKSIZE;i++){
        maskquad[i]=new float[QUAD_MASKSIZE];
        if(maskquad[i]==NULL){
             LOG(ERROR)<<"***Out of Memory Error***";
             exit(1);
        }
    }
}

Mask::Mask(char* maskfile){
    strcpy(this->maskfile,maskfile); 
     maskquad=new float*[QUAD_MASKSIZE];
    if(maskquad==NULL){
        LOG(ERROR)<<"***Out of Memory Error***";
        exit(1);
    }
    for(int i=0;i<QUAD_MASKSIZE;i++){
        maskquad[i]=new float[QUAD_MASKSIZE];
        if(maskquad[i]==NULL){
             LOG(ERROR)<<"***Out of Memory Error***";
             exit(1);
        }
    }
}

Mask::Mask(const Mask& m){
     strcpy(maskfile,m.maskfile);
     maskquad=new float*[QUAD_MASKSIZE];
     if(maskquad==NULL){
        LOG(ERROR)<<"***Out of Memory Error***";
        exit(1);
     }
     for(int i=0;i<QUAD_MASKSIZE;i++){
        maskquad[i]=new float[QUAD_MASKSIZE];
        if(maskquad[i]==NULL){
            LOG(ERROR)<<"***Out of Memory Error***";
            exit(1);
        }
     }
}
Mask::~Mask(){
    for(int i=0;i<QUAD_MASKSIZE;i++)
       delete[] maskquad[i];
    delete[] maskquad;   
}

int Mask::getMaskDet(int detid,float **mask){
    //LOG(ERROR)<<"Inside getMaskDet() - ExpMap.cpp";   
    int x1,y1,i,j,x,y;
    int posx,posy;           //position (x,y) of detector in the quadrant
    if(detid<0){
        LOG(ERROR)<<"***"<<detid<<" :Detector id cannot be negative***";
        return (EXIT_FAILURE);
    }
    posx=detid%NUM_DET_PER_QUAD_PER_ROW;
    posy=detid/NUM_DET_PER_QUAD_PER_ROW;
   
    x1=posx*(DET_MASKSIZE+NPIXELS_GAP);
    y1=posy*(DET_MASKSIZE+NPIXELS_GAP);
    
    //cout<<"x1:"<<x1<<"  y1:"<<y1<<"  posx:"<<posx<<"  posy:"<<posy;
    
    for(i=y1,y=0;y<DET_MASKSIZE;i++,y++){
        for(j=x1,x=0;x<DET_MASKSIZE;j++,x++){
            mask[y][x]=maskquad[i][j];
        }
    }
    return (EXIT_SUCCESS);
}

int Mask::getMaskQuad(int quadid){
    //cout<<"Inside getMaskQuad() - ExpMap.cpp";   
    fitsfile *fptr;
    int status=0;
    fits_open_file(&fptr,maskfile,READONLY,&status); 
    if(status){
        LOG(ERROR)<<"***Error in opening file "<<maskfile;
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    int hdutype,hdunum=0;
    fits_get_num_hdus(fptr,&hdunum,&status);  
    if(status){
        LOG(ERROR)<<"***Error getting number of HDUs in file "<<maskfile<<"***";
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    
    int qid,maxdim=2,bitpix=0,naxis=0,datatype,i;
    long fpixel[2];  fpixel[0]=fpixel[1]=1;
    long naxes[2];
    char extname[FLEN_VALUE];
    char quadname[FLEN_VALUE];
    char *keyname="EXTNAME";
    float array[QUAD_MASKSIZE*QUAD_MASKSIZE];
    //converting quadrant id to extension name for mask file 
    switch(quadid){
        case 0: strcpy(quadname,"Q0"); break;
        case 1: strcpy(quadname,"Q1"); break;
        case 2: strcpy(quadname,"Q2"); break;
        case 3: strcpy(quadname,"Q3"); break;
        default: LOG(ERROR)<<"***Invalid Quadrant ID***";
        return (EXIT_FAILURE);
    }
    //cout<<"Extension Name:"<<quadname;
    
    for(i=1;i<=hdunum;i++){
        status=0; fits_movabs_hdu(fptr,i,&hdutype,&status);
        if(status){
            LOG(ERROR)<<"***Error moving to HDU "<<i<<" in file "<<maskfile<<"***";
            return status;
        }
        status=0; fits_read_key(fptr,TSTRING,keyname,extname,NULL,&status);
        if(status!=0) continue;
        //cout<<"Extname from header:"<<extname;
        //cout<<"strcmp(extname,quadname):"<<strcmp(extname,quadname);
        if(strcmp(extname,quadname)==0 && hdutype==IMAGE_HDU){
            status=0; fits_get_img_param(fptr,maxdim,&bitpix,&naxis,naxes,&status);
            if(status){
                LOG(ERROR)<<"***Image parameters not found in file "<<maskfile<<"***";
                fits_report_error(stderr,status);
                return (EXIT_FAILURE);
            }
            if(naxes[0]!=QUAD_MASKSIZE || naxes[1]!=QUAD_MASKSIZE){
                LOG(ERROR)<<"***Dimensions for mask pattern matrix are not correct***";
                return (EXIT_FAILURE);
            }
            if(bitpix==FLOAT_IMG){
                status=0; 
                fits_read_pix(fptr,TFLOAT,fpixel,QUAD_MASKSIZE*QUAD_MASKSIZE,NULL,array,NULL,&status);
                if(status){
                    LOG(ERROR)<<"***Error reading mask data from "<<maskfile<<"***";
                    fits_report_error(stderr,status);
                    return (EXIT_FAILURE);
                }
                //cout<<"Mask pixels Reading complete";
            }
            else{
               LOG(ERROR)<<"***Mask Pattern Matrix Data Type is not Float: Must Be Float***";
               return (EXIT_FAILURE);
            }
            
            int index=0;
            for(int y=0;y<QUAD_MASKSIZE;y++){
                for(int x=0;x<QUAD_MASKSIZE;x++){
                    maskquad[y][x]=array[index++];
                    //cout<<array[index-1];
                }
            }
            break;
        }
       
    }
    fits_close_file(fptr,&status);
    if(status){
        fits_report_error(stderr,status);
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

//function to shift the detector mask array by dx and dy
//int shiftMaskDet(float mask[DET_MASKSIZE][DET_MASKSIZE],double delta_x,double delta_y){
int shiftMaskDet(float **mask,double delta_x,double delta_y){  
//cout<<"Inside shiftMaskDet() - ExpMap.cpp";   
    delta_x=delta_x*1000,delta_y=delta_y*1000;                 //converting to micrometer       
    int shiftx_pix=(delta_x)<0 ? (delta_x*(-1)/PIX_SIZ_MASK) : (delta_x/PIX_SIZ_MASK);  //shift in pixels
    int shifty_pix=(delta_y)<0 ? (delta_y*(-1)/PIX_SIZ_MASK) : (delta_y/PIX_SIZ_MASK);  //shift in pixels
                               //for one detector
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
        colstart=shiftx_pix;
        colend=DET_MASKSIZE;
        xstart=0;
    }
    else{
        colstart=0;
        colend=DET_MASKSIZE-shiftx_pix;
        xstart=shiftx_pix;
    }

    if(delta_y>0){
        rowstart=shifty_pix;
        rowend=DET_MASKSIZE;
        ystart=0;
    }
    else{
        rowstart=0;
        rowend=DET_MASKSIZE-shifty_pix;
        ystart=shifty_pix;
    }

    for(i=rowstart,y=ystart;i<rowend;i++,y++){
        for(j=colstart,x=xstart;j<colend;j++,x++){
            shiftedmask[i][j]=mask[y][x];
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
    
    for(i=0;i<DET_MASKSIZE;i++)  delete[] shiftedmask[i];
    delete[] shiftedmask;
    
    return (EXIT_SUCCESS);
}
    
//int getOpenFraction(float shiftedmask[DET_MASKSIZE][DET_MASKSIZE],
//        float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET]){


int getOpenFraction(float **shiftedmask,float **openFraction){
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

int getWeightDet(float **openFraction,float **weight){
    //cout<<"Inside getWeightDet() - ExpMap.cpp"; 
    //ofstream fout("weight.txt");
    int i,j;
    for(i=0;i<PIXELS_PER_COL_DET;i++){
        for(j=0;j<PIXELS_PER_COL_DET;j++){
            if(openFraction[i][j]==0) weight[i][j]=0;
            else  weight[i][j]=2*openFraction[i][j]-1;
            //fout<<endl<<openFraction[i][j]<<"  "<<weight[i][j];
        }
    }
    //fout.close();
    return (EXIT_SUCCESS);
}
  
void getArea(float **area,double *totarea){
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
                case 127:   l=2.28; //in mm
                break;
                default: l=2.46;  
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
                case 127:   b=2.28;     //in mm
                break;
                default:  b=2.46;
            }
            area[i][j]=l*b;
        }
    } 
    *totarea=DETECTORWIDTH*DETECTORWIDTH*TOTAL_DET_MODULES;                    //in sq mm
}

//int getWeight(float maskweight[YSIZE][XSIZE],float openf[YSIZE][XSIZE],
//        char *maskfile,float dx,float dy){
//    //cout<<"Inside getWeight() - ExpMap.cpp";   
//    int i,j,k,u,v,w;
//    Mask objmask(maskfile);
//    float detmask[DET_MASKSIZE][DET_MASKSIZE];
//    float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET];
//    float weightdet[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET];
//    unsigned char row,col; 
//    
//    for(w=0;w<NUMQUAD;w++)
//    {   
//        if(objmask.getMaskQuad(w)){             //reads mask for quadrant
//            cout<<endl<<"***Error in getting mask for quadrant "<<w<<"***"<<endl;
//            return (EXIT_FAILURE);
//        }
//        for(i=0;i<NUM_DET_PER_QUADRANT;i++){
//            if(objmask.getMaskDet(i,detmask)){
//                cout<<endl<<"***Error in getting mask for detector "<<i<<"***"<<endl;
//                return (EXIT_FAILURE);
//            }
//            if(shiftMaskDet(detmask,dx,dy)){        //shifts the mask with dx and dy
//                cout<<"***Error in shifting mask for detector by "<<dx<<","<<dy<<"***"<<endl;
//                return (EXIT_FAILURE);
//            }
//            getOpenFraction(detmask,openFraction);
//            getWeightDet(openFraction,weightdet);
//            getxy(w,i,0,&col,&row);             //w is for quad id, i is detid, 0 is pixid
//            //row=(i/NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_ROW_DET;
//            //col=(i%NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_COL_DET;
//            for(j=0,u=row;j<PIXELS_PER_COL_DET;u++,j++){
//                for(k=0,v=col;k<PIXELS_PER_ROW_DET;k++,v++){
//                    maskweight[u][v]=weightdet[j][k];
//                    openf[u][v]=openFraction[j][k];
//                    if(maskweight[u][v]==NAN || maskweight[u][v]==(-1*NAN)){
//                        cerr<<"***Mask Weight is NAN at "<<u<<","<<v<<"***";
//                        return (EXIT_FAILURE);
//                    }
//                    if(openf[u][v]==NAN){
//                        cerr<<"***Open Fraction is NAN at "<<u<<","<<v<<"***";
//                        return (EXIT_FAILURE);
//                    }
//                    //cout<<endl<<u<<","<<v<<" - "<<maskweight[u][v]<<" "<<openf[u][v];
//                }
//            }
//        }
//    }
//    return (EXIT_SUCCESS);
//}

int cztgenshadow(double thetax,double thetay,double sourceflux,char *maskfile,
        float shadow[YSIZE][XSIZE]){
    //cout<<"Inside cztgenshadow() - ExpMap.cpp";   
    double dx,dy;
    dx=MASKHEIGHT*tan(thetax);
    dy=MASKHEIGHT*tan(thetay);
    LOG(INFO)<<"dx:"<<dx<<"   dy:"<<dy;
    if(dx>DETECTORWIDTH || dy>DETECTORWIDTH){
        LOG(ERROR)<<"***Shift greater than the limit***";
        return (EXIT_FAILURE);
    }
    float **weight,**openf;
    float **ai;
    
    weight=new float*[YSIZE];
    openf=new float*[YSIZE];
    ai=new float*[YSIZE];
   
    if(weight==NULL || openf==NULL || ai==NULL){
        LOG(ERROR)<<"***Out of Memory***";
        return (EXIT_FAILURE);
    }
    for(int i=0;i<YSIZE;i++){
            weight[i]=new float[XSIZE];
            openf[i]=new float [XSIZE];
            ai[i]=new float[XSIZE];
            if(weight[i]==NULL || openf[i]==NULL || ai[i]==NULL){
                LOG(ERROR)<<"***Out of Memory***";
                return (EXIT_FAILURE);
            }
     }
    
    
     if(getWeight(weight,openf,maskfile,dx,dy)){
        LOG(ERROR)<<"***Error in getWeight***";
        return (EXIT_FAILURE);
    }
    
    double totarea=0;
    getArea(ai,&totarea);
    if(totarea==0){
        LOG(ERROR)<<"***Total Area="<<totarea<<" Divide by zero error***";
        return (EXIT_FAILURE);
    }
    double C=sourceflux/totarea;
    double sourceangle=atan(sqrt(tan(thetax)*tan(thetax)+tan(thetay)*tan(thetay)));
    
    //cout<<"shadow";
    for(int i=0;i<YSIZE;i++){
        for(int j=0;j<XSIZE;j++){
            shadow[i][j]=C*cos(sourceangle)*openf[i][j]*ai[i][j];
            //cout<<"|"<<shadow[i][j];
        }
    }
    
    for(int i=0;i<YSIZE;i++){
        delete[] weight[i],openf[i],ai[i];
    }
    delete[] weight,openf,ai;
    
    return (EXIT_SUCCESS);
}
                                   
int getExpMap(double theta_x,double theta_y,char *maskfile,float emap[XSIZE*YSIZE]){
    LOG(INFO)<<"Creating Exposure Map.........";
    float emap_2d[YSIZE][XSIZE];
    double dx,dy;
    dx=MASKHEIGHT*tan(theta_x);
    dy=MASKHEIGHT*tan(theta_y);
    LOG(INFO)<<"dx:"<<dx<<"   dy:"<<dy;
    if(dx>DETECTORWIDTH || dy>DETECTORWIDTH){
        LOG(ERROR)<<"***Shift greater than the limit***";
        return (EXIT_FAILURE);
    }
    //------
    int i,j,k,u,v,w;
    Mask objmask(maskfile);
    //float detmask[DET_MASKSIZE][DET_MASKSIZE];
    //float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET];
    
    float **detmask,**openFraction;
    openFraction=new float *[PIXELS_PER_COL_DET];
    detmask=new float *[DET_MASKSIZE];
    if(openFraction==NULL || detmask==NULL){
         LOG(ERROR)<<"Out of memory";
         return (EXIT_FAILURE);
    }
    for(int i=0;i<PIXELS_PER_COL_DET;i++){
        openFraction[i]=new float[PIXELS_PER_ROW_DET];
        if(openFraction[i]==NULL){
            LOG(ERROR)<<"Out of memory";
            return (EXIT_FAILURE);
        }
    }
    for(int i=0;i<DET_MASKSIZE;i++){
        detmask[i]=new float[DET_MASKSIZE];
        if(detmask[i]==NULL){
            LOG(ERROR)<<"Out of memory";
            return (EXIT_FAILURE);
        }
    }
    
    
    
    unsigned char row,col; 
    for(w=0;w<NUMQUAD;w++)
    {   LOG(ERROR)<<"Quadrant Number-"<<w;
        if(objmask.getMaskQuad(w)){             //reads mask for quadrant
            LOG(ERROR)<<"***Error in getting mask for quadrant "<<w<<"***";
            return (EXIT_FAILURE);
        }
        for(i=0;i<NUM_DET_PER_QUADRANT;i++){
            objmask.getMaskDet(i,detmask);
            if(shiftMaskDet(detmask,dx,dy)){        //shifts the mask with dx and dy
                LOG(ERROR)<<"***Error in shifting mask for detector by "<<dx<<","<<dy;
                return (EXIT_FAILURE);
            }
            getOpenFraction(detmask,openFraction);
            getxy(w,i,0,&col,&row);
            //row=(i/NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_ROW_DET;
            //col=(i%NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_COL_DET;
            for(j=0,u=row;j<PIXELS_PER_COL_DET;u++,j++){
                for(k=0,v=col;k<PIXELS_PER_ROW_DET;k++,v++){
                    //fout<<""<<setw(6)<<u<<setw(6)<<v<<setw(6)<<j<<setw(6)<<k;
                    emap_2d[u][v]=openFraction[j][k];
                    if(emap_2d[u][v]==NAN){
                        LOG(ERROR)<<"***Open Fraction is NAN at "<<u<<","<<v;
                        return (EXIT_FAILURE);
                    }
                    //cout<<endl<<u<<","<<v<<" - "<<maskweight[u][v]<<" "<<openf[u][v];
                }
            }
        }
    }
    
    //writeImg("emap_2d",emap_2d,XSIZE,YSIZE);
    
    k=0;
    for(i=0;i<YSIZE;i++)
        for(j=0;j<XSIZE;j++)
            emap[k++]=emap_2d[i][j];
     for(int i=0;i<PIXELS_PER_COL_DET;i++){
        delete[] openFraction[i];
    }
    delete[] openFraction;
    
    for(int i=0;i<DET_MASKSIZE;i++)
        delete[] detmask[i];
    delete[] detmask;
    
   
    LOG(INFO)<<"Created Exposure Map";
    return (EXIT_SUCCESS);
}

//int writeImg(char *file,float img[XSIZE][YSIZE],int m,int n){
//    float data[m*n];
//    int index=0;
//    for(int i=0;i<m;i++)
//        for(int j=0;j<n;j++)
//            data[index++]=img[i][j];
//    int status=0;
//    fitsfile *fptr;
//    fits_create_file(&fptr,file,&status);
//    int bitpix=FLOAT_IMG;
//    int naxis=2;
//    long naxes[2];
//    naxes[0]=m; naxes[1]=n;
//    fits_create_img(fptr,bitpix,naxis,naxes,&status);
//    long fpixel[2];
//    fpixel[0]=fpixel[1]=1;
//    fits_write_pix(fptr,TFLOAT,fpixel,m*n,data,&status);
//    fits_close_file(fptr,&status);
//    return 0;        
//}


int getWeight(float **maskweight,float **openf,char *maskfile,float dx,float dy){
    //cout<<"Inside getWeight() - ExpMap.cpp";   
    int i,j,k,u,v,w;
    Mask objmask(maskfile);
//    float detmask[DET_MASKSIZE][DET_MASKSIZE];
//    float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET];
//    float weightdet[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET];
    
    float **detmask,**openFraction,**weightdet;
    openFraction=new float *[PIXELS_PER_COL_DET];
    weightdet=new float *[PIXELS_PER_COL_DET];
    detmask=new float *[DET_MASKSIZE];
    if(openFraction==NULL || weightdet==NULL || detmask==NULL){
         LOG(ERROR)<<"Out of memory";
         return (EXIT_FAILURE);
    }
    for(int i=0;i<PIXELS_PER_COL_DET;i++){
        openFraction[i]=new float[PIXELS_PER_ROW_DET];
        weightdet[i]=new float[PIXELS_PER_ROW_DET];
        if(weightdet[i]==NULL || openFraction[i]==NULL){
            LOG(ERROR)<<"Out of memory";
            return (EXIT_FAILURE);
        }
    }
    for(int i=0;i<DET_MASKSIZE;i++){
        detmask[i]=new float[DET_MASKSIZE];
        if(detmask[i]==NULL){
            LOG(ERROR)<<"Out of memory";
            return (EXIT_FAILURE);
        }
    }
        
    
    unsigned char row,col; 
    
    for(w=0;w<NUMQUAD;w++)
    {   
        if(objmask.getMaskQuad(w)){             //reads mask for quadrant
            LOG(ERROR)<<"***Error in getting mask for quadrant "<<w<<"***";
            return (EXIT_FAILURE);
        }
        for(i=0;i<NUM_DET_PER_QUADRANT;i++){
            if(objmask.getMaskDet(i,detmask)){
                LOG(ERROR)<<"***Error in getting mask for detector "<<i<<"***";
                return (EXIT_FAILURE);
            }
            if(shiftMaskDet(detmask,dx,dy)){        //shifts the mask with dx and dy
                LOG(ERROR)<<"***Error in shifting mask for detector by "<<dx<<","<<dy<<"***";
                return (EXIT_FAILURE);
            }
            getOpenFraction(detmask,openFraction);
            getWeightDet(openFraction,weightdet);
            getxy(w,i,0,&col,&row);             //w is for quad id, i is detid, 0 is pixid
            //row=(i/NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_ROW_DET;
            //col=(i%NUM_DET_PER_QUAD_PER_ROW)*PIXELS_PER_COL_DET;
            for(j=0,u=row;j<PIXELS_PER_COL_DET;u++,j++){
                for(k=0,v=col;k<PIXELS_PER_ROW_DET;k++,v++){
                    maskweight[u][v]=weightdet[j][k];
                    openf[u][v]=openFraction[j][k];
                    if(maskweight[u][v]==NAN || maskweight[u][v]==(-1*NAN)){
                        LOG(ERROR)<<"***Mask Weight is NAN at "<<u<<","<<v<<"***";
                        return (EXIT_FAILURE);
                    }
                    if(openf[u][v]==NAN){
                        LOG(ERROR)<<"***Open Fraction is NAN at "<<u<<","<<v<<"***";
                        return (EXIT_FAILURE);
                    }
                    //cout<<u<<","<<v<<" - "<<maskweight[u][v]<<" "<<openf[u][v];
                }
            }
        }
    }
    
    for(int i=0;i<PIXELS_PER_COL_DET;i++){
        delete[] openFraction[i];
        delete[] weightdet[i];
    }
    delete[] openFraction, weightdet;
    
    for(int i=0;i<DET_MASKSIZE;i++)
        delete[] detmask[i];
    delete[] detmask;
    
    return (EXIT_SUCCESS);
}
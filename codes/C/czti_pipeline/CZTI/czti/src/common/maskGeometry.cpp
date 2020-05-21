#include "maskGeometry.h"

TanMaskGeometry::TanMaskGeometry(){
    
}

int TanMaskGeometry::init_full_finer_mask() {
    int status=0;
    int i=0;
    
    //initializing vectors xmask0 & ymask0
    xmask0.clear();
    ymask0.clear();
    for(i=0; i<NY_MASKPIXELS; i++){
        xmask0.push_back(0.0);
        ymask0.push_back(0.0);
    }
    
    for(i=0; i<NY_MASKPIXELS; i++){
        xmask0[i]=i*SUB_PIXEL_WIDTH;
        ymask0[i]=i*SUB_PIXEL_WIDTH;
    }
    
    return EXIT_SUCCESS;
}



int TanMaskGeometry::read_compressed_mask(string compMaskFilename) {
    int status=0;
    int i=0, j=0, k=0;
    fitsfile *fcompmask;
    long fpixel[2]; fpixel[0]=1; fpixel[1]=1;
    long npixels = COMPROWS*COMPCOLS;
    char hduname[25];
    unsigned int* compMaskArray; //Compressed mask for a quadrant
    unsigned char* unFieldArray; // uncompressed field values (25)
    unsigned char* unMaskArray;  // Uncompressed Quadrant mask (finer)
    vector <unsigned char> vecTempUnMask; 
    
    
    //Initializing variables
    compMaskArray = new unsigned int [npixels]; //Compressed mask array for a single quadrant
    unFieldArray = new unsigned char [COMP_FIELD_SIZE]; //uncompressed field consisting of 25 values
    unMaskArray = new unsigned char [QUAD_MASKSIZE*QUAD_MASKSIZE];
    for(i=0; i<npixels; i++){
        compMaskArray[i]=0;
    }
    for(i=0; i<COMP_FIELD_SIZE; i++){
        unFieldArray[i]=0;
    }
    for(i=0; i<QUAD_MASKSIZE*QUAD_MASKSIZE; i++){
        unMaskArray[i]=0; 
    }
    
    vecTempUnMask.resize(QUAD_MASKSIZE, 0);
    
    q0UnMask.clear();
    q1UnMask.clear();
    q2UnMask.clear();
    q3UnMask.clear();
    q0UnMask.resize(QUAD_MASKSIZE,vecTempUnMask);
    q1UnMask.resize(QUAD_MASKSIZE,vecTempUnMask);
    q2UnMask.resize(QUAD_MASKSIZE,vecTempUnMask);
    q3UnMask.resize(QUAD_MASKSIZE,vecTempUnMask);
    //Variables initialized

    DLOG(INFO) << "Extracting fine mask pattern from compressed mask file: " << compMaskFilename;
    
    //Opening and reading compressed mask file.
    fits_open_file(&fcompmask, (char *) compMaskFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening compressed mask file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    for(i=0; i<NUMQUAD; i++){
        quadToHDU(i, hduname);
        fits_movnam_hdu(fcompmask, IMAGE_HDU, hduname, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << hduname << " in compressed mask file: " << compMaskFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        fits_read_pix(fcompmask, TUINT, fpixel, npixels, NULL, compMaskArray, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading compressed mask for " << hduname << " from compressed mask file: " << compMaskFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        for(j=0; j<npixels; j++){
            //DLOG(INFO) << j;
            uncompress(compMaskArray, unFieldArray, j);
            for(k=0; k<COMP_FIELD_SIZE; k++){
                unMaskArray[j*COMP_FIELD_SIZE+k] = unFieldArray[k];
                //cout  <<  (int) unMaskArray[j*COMP_FIELD_SIZE + k] << "  ";
                
            }
            //cout << endl;
        }
        DLOG(INFO)<< "Quadrant " << i << " mask has been uncompressed.";
        
        //Storing data in 2D vectors
        for(j=0; j<QUAD_MASKSIZE; j++){
            for(k=0; k<QUAD_MASKSIZE; k++) {
                if (i == 0) {q0UnMask[j][k] = unMaskArray[j * QUAD_MASKSIZE + k];}
                else if (i == 1) {q1UnMask[j][k] = unMaskArray[j * QUAD_MASKSIZE + k];}
                else if (i == 2) {q2UnMask[j][k] = unMaskArray[j * QUAD_MASKSIZE + k];}
                else if (i == 3) {q3UnMask[j][k] = unMaskArray[j * QUAD_MASKSIZE + k];}
            }
        }    
    }//END LOOP ON UNCOMPRESSED MASK HDUS
    
    fits_close_file(fcompmask, &status);
    if (status) {
        LOG(ERROR) << "Error in closing compressed mask file: " << compMaskFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
     
    delete[] compMaskArray, unFieldArray, unMaskArray;
    return EXIT_SUCCESS;
}

int TanMaskGeometry::uncompress(unsigned int *maskpattern,unsigned char *uncompressed_field, int i)
{
        int k=0,val=0, temp_val=0;
        //DLOG(INFO) << "Mask Pattern [" << i << "] : " << maskpattern[i] << endl;
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

int TanMaskGeometry::read_64x64_mask(string maskFilename){
    int status=0;
    int i=0, k=0; //counter variables
    int j=0;
    int anynull;
    long naxes[2], nbuffer, npixels=XPIX_QUAD*YPIX_QUAD;
    long fpixel[2];
    float nullval=0;
    fitsfile *fmask;
    char hduname[25];
    vector <unsigned char> vecTempMask;
    float* quadMaskArray;
    
    fpixel[0]=1;
    fpixel[1]=1;
    //clearing class variables
    q0Mask.clear();
    q1Mask.clear();
    q2Mask.clear();
    q3Mask.clear();
    fullMask.clear();
    
    //Initializing variables
    quadMaskArray = new float [XPIX_QUAD*YPIX_QUAD];
    for(i=0; i<XPIX_QUAD*YPIX_QUAD; i++){
            quadMaskArray[i]=0.0;
    }
    
    vecTempMask.clear();
    for(i=0; i<XPIX_QUAD; i++){
        vecTempMask.push_back(0);
    }
    for(i=0; i<YPIX_QUAD; i++){
        q0Mask.push_back(vecTempMask);
        q1Mask.push_back(vecTempMask);
        q2Mask.push_back(vecTempMask);
        q3Mask.push_back(vecTempMask);
    }
    //Opening fits file to read mask 
    fits_open_file(&fmask, (char*) maskFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening CALDB 64x64 mask file " << maskFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    for(i=0; i<NUMQUAD; i++){
        
        DLOG(INFO) << "Reading mask for Quadrant " << i << "...";
        quadToHDU(i, hduname);
        fits_movnam_hdu(fmask, IMAGE_HDU, hduname, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << hduname << " in mask file: " << maskFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_read_pix(fmask, TFLOAT, fpixel, npixels, NULL, quadMaskArray, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading mask for " << hduname << " from mask file: " << maskFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        
        for(j=0; j<YPIX_QUAD; j++){
            for(k=0; k<XPIX_QUAD; k++){
                if(i==0){q0Mask[j][k]=(unsigned char) quadMaskArray[j*YPIX_QUAD + k];}
                else if(i==1){q1Mask[j][k]=(unsigned char) quadMaskArray[j*YPIX_QUAD + k];}
                else if(i==2){q2Mask[j][k]=(unsigned char) quadMaskArray[j*YPIX_QUAD + k];}
                else if(i==3){q3Mask[j][k]=(unsigned char) quadMaskArray[j*YPIX_QUAD + k];}
            }
        }

    }

    fits_close_file(fmask, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file: " << maskFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
   
    delete[] quadMaskArray; 
    
    return EXIT_SUCCESS;
}

int TanMaskGeometry::get_full_finer_mask(){
    int status=0;
    int i,j=0;
    vector <unsigned char> vecTempMask;
    
    //Checking whether full mask can be generated or not
    if(q0UnMask.size()==0 || q1UnMask.size()==0 ||q2UnMask.size()==0 || q3UnMask.size()==0) {
            LOG(ERROR) << "Full 64x64 mask can't be created as one of the quadrant mask has NULL size.";
    }

    //Initializing class variables
    vecTempMask.clear();
    fullUnMask.clear();
    for (i = 0; i < (QUAD_MASKSIZE*2 + PIX_DIST_QUADRANTS); i++) {
        vecTempMask.push_back(1);
    }
    for (i = 0; i < (QUAD_MASKSIZE*2 + PIX_DIST_QUADRANTS); i++) {
        fullUnMask.push_back(vecTempMask);
    }
    DLOG(INFO) << "Size of full mask: " << vecTempMask.size() << "x" << fullUnMask.size();

    if (rearrange_finer_mask_quads(&fullUnMask, q0UnMask, 0)) {
        LOG(ERROR) << "Error in properly placing q0mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_finer_mask_quads(&fullUnMask, q1UnMask, 1)) {
        LOG(ERROR) << "Error in properly placing q1mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_finer_mask_quads(&fullUnMask, q2UnMask, 2)) {
        LOG(ERROR) << "Error in properly placing q2mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_finer_mask_quads(&fullUnMask, q3UnMask, 3)) {
        LOG(ERROR) << "Error in properly placing q3mask into full mask";
        return EXIT_FAILURE;
    }    
    DLOG(INFO) << "Full finer mask constructed from individual quadrant uncompressed masks";
    return EXIT_SUCCESS;
    
}
int TanMaskGeometry::get_full64x64_mask(){
    int status=0;
    int i,j=0;
    vector <unsigned char> vecTempMask;
    
    //Checking whether full mask can be generated or not
    if(q0Mask.size()==0 || q1Mask.size()==0 || q2Mask.size()==0 || q3Mask.size()==0){
        LOG(ERROR) << "Full 64x64 mask can't be created as one of the quadrant mask has NULL size.";
    }
    
    //Initializing class variables
    vecTempMask.clear();
    fullMask.clear();
    for(i=0; i<XSIZE; i++){
        vecTempMask.push_back(0);
    }
    for(i=0; i<YSIZE; i++){
        fullMask.push_back(vecTempMask);
    }

    if (rearrange_quads(&fullMask, q0Mask, 0)) {
        LOG(ERROR) << "Error in properly placing q0mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_quads(&fullMask, q1Mask, 1)) {
        LOG(ERROR) << "Error in properly placing q1mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_quads(&fullMask, q2Mask, 2)) {
        LOG(ERROR) << "Error in properly placing q2mask into full mask";
        return EXIT_FAILURE;
    }
    if (rearrange_quads(&fullMask, q3Mask, 3)) {
        LOG(ERROR) << "Error in properly placing q3mask into full mask";
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

TanMaskWeighting::TanMaskWeighting(double thetaX, double thetaY, float flux) {
    int i=0, j=0;
    this->thetaX = thetaX;
    this->thetaY = thetaY;
    vector <float> tempOpenFraction;
    //Initializing vectors
    for(i=0; i<PIXELS_PER_ROW_DET; i++){
        tempOpenFraction.push_back(0.0);
    }
    for(i=0; i<PIXELS_PER_COL_DET; i++){
        openFractionDet.push_back(tempOpenFraction);
    }
    
    //computing mask pattern shift in X and Y direction in cm.
    if(compute_shift(thetaX, thetaY, this->dX, this->dY)){
        LOG(ERROR) << "Error in computing shift.";
    }
}

int TanMaskWeighting::calculate_open_fraction_det(vector<vector<unsigned char> > detMask) {
    int ipixx=0, ipixy=0, x=0, y=0;
    int xstart, xend, ystart, yend;
    int wpixSmall, wpixBig; //width of small and big pixels
    int closecount=0; //number of closed pixels
    float totpix=0; //number of 
    
    wpixSmall = (SMALL_PIXEL_WIDTH/PIX_SIZ_MASK) - 1;
    wpixBig = (BIG_PIXEL_WIDTH/PIX_SIZ_MASK);
    DLOG(INFO) << "width_pix_small: " << wpixSmall;
    DLOG(INFO) << "width_pix_big: " << wpixBig;
    
    for(ipixy=0; ipixy<PIXELS_PER_COL_DET; ipixy++){
        for(ipixx=0; ipixx<PIXELS_PER_ROW_DET; ipixx++){
            closecount=0;
            if(ipixx==0) {xstart=0; xend=wpixSmall;}
            else if(ipixx==15) {xstart=(wpixSmall+14*wpixBig); xend=xstart+wpixSmall;}
            else {xstart=(wpixSmall+(ipixx-1)*wpixBig); xend=(xstart+wpixBig);}
            
            if(ipixy==0) {ystart=0; yend=wpixSmall;}
            else if(ipixy==15) {ystart=(wpixSmall+14*wpixBig); yend=ystart+wpixSmall;}
            else {ystart=(wpixSmall+(ipixy-1)*wpixBig); yend=(ystart+wpixBig);}

            if ((ipixx==0 && ipixy==0) || (ipixx==15 && ipixy==0) || (ipixx==0 && ipixy==15) || (ipixx==15 && ipixy==15)) totpix = wpixSmall*wpixSmall;
            else if ((ipixx==0 && 0 < ipixy < 15) || (ipixx==15 && 0 < ipixy < 15) || (ipixy==0 && 0 < ipixx < 15) || (ipixy==15 && 0 < ipixx < 15)) totpix = wpixSmall * wpixBig;
            else totpix = wpixBig*wpixBig;

            for (y = ystart; y < yend; y++) {
                for (x = xstart; x < xend; x++) {
                    if (detMask[y][x] == CLOSE)
                        closecount++;
                }
            }
            openFractionDet[ipixy][ipixx] = 1.0 - ((float) closecount / (float) totpix); 
            //DLOG(INFO) << "detector (i,j): (" << ipixx << "," << ipixy << ")";
            //DLOG(INFO) << "closecount: " << closecount;
            //DLOG(INFO) << "totpix: " << totpix;
            //DLOG(INFO) << "openfraction: " << openFractionDet[ipixy][ipixx];
        }
    }
    return EXIT_SUCCESS;
}

int TanMaskWeighting::compute_open_fraction_quad(int quadNo, TanMaskGeometry tanMask){
    int status=0;
    int i=0;
    int detid=0;
    vector < vector <unsigned char> > shiftedDetMask;
    vector <float> tempOpenFraction;
    
    
    //Initializing vectors
    for(i=0; i<XPIX_QUAD; i++){
        tempOpenFraction.push_back(0.0);
    }
    if(quadNo==0){
        openFractionQ0.clear();
        for(i=0; i<YPIX_QUAD; i++){
            openFractionQ0.push_back(tempOpenFraction);
        }
    }
    else if(quadNo==1){
        openFractionQ1.clear();
        for(i=0; i<YPIX_QUAD; i++){
            openFractionQ1.push_back(tempOpenFraction);
        }
    }
    else if(quadNo==2){
        openFractionQ2.clear();
        for(i=0; i<YPIX_QUAD; i++){
            openFractionQ2.push_back(tempOpenFraction);
        }
    }
    else if(quadNo==3){
        openFractionQ3.clear();
        for(i=0; i<YPIX_QUAD; i++){
            openFractionQ3.push_back(tempOpenFraction);
        }
    }
    else {
        LOG(ERROR) << "Error in Quadrant Number: " << quadNo << ". Quadrant number should lie in the range 0-3.";
        return (EXIT_FAILURE);
    }
    //vectors initialized
    
    DLOG(INFO) << "Computing openfraction map for quadrant " << quadNo;
    for(detid=0; detid<NO_DET_PER_QUAD; detid++){
        shiftedDetMask = get_shifted_detector_mask(quadNo, detid, tanMask, dX, dY, status);
        if(status){
            LOG(ERROR) << "Error in getting shifted detector mask for detector id " << detid << " of quadrant " << quadNo;
            return (EXIT_FAILURE);
        }
        
        //calculating openfraction of a particular detector by using shifted detector mask
        if(calculate_open_fraction_det(shiftedDetMask)){
            LOG(ERROR) << "Error in calculating open fraction for detector id " << detid << " of quadrant " << quadNo;
            return (EXIT_FAILURE);
        }
        DLOG(INFO) << "[DEBUG]";
        //Placing detector open fraction in relevant quadrant open fraction.
        if(quadNo==0) { 
            if(rearrange_detectors(&openFractionQ0, openFractionDet, detid)){
                LOG(ERROR) << "Error in generation of Quadrant " << quadNo << " open fraction.";
                return (EXIT_FAILURE);
            }
        }
        else if(quadNo==1) { 
            if(rearrange_detectors(&openFractionQ1, openFractionDet, detid)){
                LOG(ERROR) << "Error in generation of Quadrant " << quadNo << " open fraction.";
                return (EXIT_FAILURE);
            }
        }
        else if(quadNo==2) { 
            if(rearrange_detectors(&openFractionQ2, openFractionDet, detid)){
                LOG(ERROR) << "Error in generation of Quadrant " << quadNo << " open fraction.";
                return (EXIT_FAILURE);
            }
        }
        else if(quadNo==3) { 
            if(rearrange_detectors(&openFractionQ0, openFractionDet, detid)){
                LOG(ERROR) << "Error in generation of Quadrant " << quadNo << " open fraction.";
                return (EXIT_FAILURE);
            }
        }
    }
    
    return (EXIT_SUCCESS);
}


int TanMaskWeighting::compute_open_fraction(TanMaskGeometry tanMask){
    int status=0;
    int i=0;
    vector <float> tempOpenFraction;
    
    // Initializing vectors
    for(i=0; i<XSIZE; i++){
        tempOpenFraction.push_back(0.0);
    }
    for(i=0; i<YSIZE; i++){
        openFractionFull.push_back(tempOpenFraction);
    }
    // Vectors initialized
    
    // Computing open fraction for all quadrants
    for(i=0; i<NUMQUAD; i++){
        if (compute_open_fraction_quad(i, tanMask)){
            LOG(ERROR) << "Error in computing open fraction for Quadrant " << i;
            return (EXIT_FAILURE);
        }
        
    }
    
    if(rearrange_quads(&openFractionFull, openFractionQ0, 0)){
        LOG(ERROR)<<"Error in arranging open fraction quadrants to get full image for Quadrant 0";
        return (EXIT_FAILURE);
    }
    if(rearrange_quads(&openFractionFull, openFractionQ1, 1)){
        LOG(ERROR)<<"Error in arranging open fraction quadrants to get full image for Quadrant 1";
        return (EXIT_FAILURE);
    }
    if(rearrange_quads(&openFractionFull, openFractionQ2, 2)){
        LOG(ERROR)<<"Error in arranging open fraction quadrants to get full image for Quadrant 2";
        return (EXIT_FAILURE);
    }
    if(rearrange_quads(&openFractionFull, openFractionQ3, 3)){
        LOG(ERROR)<<"Error in arranging open fraction quadrants to get full image for Quadrant 3";
        return (EXIT_FAILURE);
    }
    
    
    return EXIT_SUCCESS;
}  

TransMask::TransMask(double thetaX, double thetaY) {
    this->thetaX = thetaX;
    this->thetaY = thetaY;
}

int TransMask::compute_mask_openfraction(CztDetectorGeometry detgeom, int oversampling, 
        float energy){
    int status=0;
    int nx=0, ny=0; //oversampling in x and y direction
    float x0=0.0, y0=0.0;
    float x=0.0, y=0.0; //actual coordinates of sub-pixels
    float offset=0.0, tt=0.0;
    float deltax=0.0, deltay =0.0; //edge length and breadth of each sub-pixel.
    long idetx=0, idety=0, ipixx=0, ipixy=0, isx=0, isy=0;
    vector <float> xdet0, ydet0, pixx0, pixy0, pixxw, pixyw;
    float xmoff=0.0, ymoff=0.0;
    int pixstate=0, nclose_pixels=0; //pixstate: state of pixel (open/close)
                                     //nclose_pixels: number of close pixels.
    float openfraction;
    float oftrans; //open fraction after considering transmission from close pixels
    float tau=0.0; //optical depth of the mask plate in the given energy range
    double trad=0.0; //source angle with camera normal in radians
    double tradx=0.0, trady=0.0; //thetax and thetay in radians
    float abscoTan;
    float xpix0, ypix0; //bottom-left sub-pixel coordinate of each pixel
    float xpixw, ypixw; //length & breadth of a pixel
    
    //Resetting oftable
    oftable.reset();
    
    //limit angle to avoid zero divide
    if(abs(90.0-abs(thetaX)) < 1.0){
        thetaX=89.50;
    }
    if(abs(90.0-abs(thetaY)) < 1.0){
        thetaY=89.50;
    }

    //Avoid degeneracy of surface identity at azimuth diagonals by applying a small offset to thetax
    offset = 0.01;
    if (thetaX == thetaY) {
        if (thetaX != 0.0) {
            tt = abs(thetaX) - abs(offset);
            if (thetaX < 0.0) {
                thetaX = -abs(tt);
            } else {
                thetaX = abs(tt);
            }

        }
    }
    
    //converting theta from degrees into radians
    tradx = thetaX * M_PI/180;
    trady = thetaY * M_PI/180;
    trad = atan(sqrt(tan(tradx)*tan(tradx) + tan(trady)*tan(trady)));

    //calculating effective absorption coefficient of tantalum
    abscoTan = absco_Ta(energy);
    if(status){
        LOG(ERROR) << "Error in calculating absorption coefficient of Tantalum "
                << "for energy " << energy;
        return EXIT_FAILURE;
    }
    
    //optical depth of the mask plate
    tau = abscoTan*TANMASK_THICKNESS/cos(trad);
            
    //computing shift
    if(compute_shift(thetaX, thetaY, dx, dy)){
        LOG(ERROR)<<"Error in computing shift";
        return (EXIT_FAILURE);
    }

    

    xdet0=detgeom.get_xdet0();
    ydet0=detgeom.get_ydet0();
    detgeom.display_detector_coordinates();
    
    nx=oversampling;
    ny=oversampling;
    //loop over detectors 8x8 array
    for(idety=0; idety<8; idety++){
        y0=ydet0[idety];
        for(idetx=0; idetx<8; idetx++){
            x0=xdet0[idetx];
            detgeom.define_pixel_coordinates(x0, y0);
            pixx0 = detgeom.get_pixx0();
            pixy0 = detgeom.get_pixy0();
            pixxw = detgeom.get_pixxw();
            pixyw = detgeom.get_pixyw();
            
            for(ipixy=0; ipixy<PIXELS_PER_COL_DET; ipixy++){
                for(ipixx=0; ipixx<PIXELS_PER_ROW_DET; ipixx++){
                    //bottom-left coordinates of each pixel
                    xpix0=pixx0[ipixx];
                    ypix0=pixy0[ipixy];
                    //edge length and breadth of each pixel
                    xpixw=pixxw[ipixx];
                    ypixw=pixyw[ipixy];
                    //edge length and breadth of each sub-pixel
                    deltax=xpixw/(1.0*nx);
                    deltay=ypixw/(1.0*ny);
                    for(isy=0; isy<ny; isy++){
                        y=ypix0+(1.0*isy + 0.5)*deltay;
                        for(isx=0; isx<nx; isx++){
                            x=xpix0+(1.0*isx + 0.5)*deltax;
                            //Here (x,y) is the centre of the current area bin.
                            xmoff=x+dx;
                            ymoff=y+dy;
                            //Here (xmoff,ymoff) is the projection of current area bin on mask.
                            transmask(xmoff, ymoff, pixstate);
                            if(pixstate==CLOSE){
                                nclose_pixels++;
                            }
                        }//END LOOP OVER X SUB-PIXEL
                    }//END LOOP OVER Y SUB-PIXEL
                    //open fraction with transmission from tantalum mask
                    openfraction = (nx*ny - nclose_pixels)/(nx*ny * 1.0);
                    
                    //open fraction taking into account tantalum mask transmission
                    oftrans = (openfraction) + (1-openfraction)*exp(-tau); 
                    oftable.pushback_openfraction(idetx, idety, ipixx, ipixy, 
                            openfraction, oftrans);
                    nclose_pixels=0;
                }//END LOOP OVER X PIXEL
            }//END LOOP OVER Y PIXEL
        }//END LOOP OVER X DETECTOR
    }// END LOOP OVER Y DETECTOR
    return (EXIT_SUCCESS);
}

int TransMask::transmask(float x, float y, int &pixstate){
    int ix, iy; //integral position of x and y on finer mask
    
    ix = (x-XMIN_MASK)/SUB_PIXEL_WIDTH;
    iy = (y-YMIN_MASK)/SUB_PIXEL_WIDTH;
    //DLOG(INFO) << "INTEGRAL POSITION ON MASK :" << ix << ", " << iy; 
    if(ix>=0 && ix<NX_MASKPIXELS && iy>=0 && iy<NY_MASKPIXELS){
        pixstate = (int) fullUnMask[iy][ix];
    }
    else{
        pixstate = OPEN; //if x,y outside mask then no tantalum blocking hence 1
    }
    return (EXIT_SUCCESS);
}

int TransMask::generate_full_openfraction_image() {
    int status=0;
    long i=0,ix=0, iy=0;
    long nrows=oftable.get_nrows(); //number of rows in mask open fraction table.
    vector <float> tempOpenFraction;
    //Initializing vectors
    openFractionFull.clear();
    for(ix=0; ix<XSIZE; ix++){
        tempOpenFraction.push_back(0.0);
    }
    for(iy=0; iy<YSIZE; iy++){
        openFractionFull.push_back(tempOpenFraction);
    }
    //Vectors initialized
    
    for(i=0; i<nrows; i++){
        //Calculating pixel co-ordinates for each open fraction table record.
        ix = oftable.detX[i]*PIXELS_PER_ROW_DET + oftable.pixX[i];
        iy = oftable.detY[i]*PIXELS_PER_COL_DET + oftable.pixY[i];
        
        //Full Open fraction 2D vector 
        openFractionFull[iy][ix]=oftable.openfrac[i];
    }
    
    return EXIT_SUCCESS;
}
int TransMask::generate_full_oftrans_image() {
    int status=0;
    long i=0,ix=0, iy=0;
    long nrows=oftable.get_nrows(); //number of rows in mask open fraction table.
    vector <float> tempOFtrans;
    //Initializing vectors
    oftransFull.clear();
    for(ix=0; ix<XSIZE; ix++){
        tempOFtrans.push_back(0.0);
    }
    for(iy=0; iy<YSIZE; iy++){
        oftransFull.push_back(tempOFtrans);
    }
    //Vectors initialized
    
    for(i=0; i<nrows; i++){
        //Calculating pixel co-ordinates for each open fraction table record.
        ix = oftable.detX[i]*PIXELS_PER_ROW_DET + oftable.pixX[i];
        iy = oftable.detY[i]*PIXELS_PER_COL_DET + oftable.pixY[i];
        
        //Full Open fraction 2D vector 
        oftransFull[iy][ix]=oftable.oftrans[i];
    }
    
    return EXIT_SUCCESS;
}


//MaskOpenFractionTable struct
int MaskOpenFractionTable::pushback_openfraction(int detx, int dety, int pixx, int pixy, 
        float openfraction, float oftransmission){
    detX.push_back(detx);
    detY.push_back(dety);
    pixX.push_back(pixx);
    pixY.push_back(pixy);
    openfrac.push_back(openfraction);
    oftrans.push_back(oftransmission);
    
    return EXIT_SUCCESS;
}

int MaskOpenFractionTable::reset() {
    detX.clear();
    detY.clear();
    pixX.clear();
    pixY.clear();
    openfrac.clear();
    oftrans.clear();

    return EXIT_SUCCESS;
}


//Friend function of TanMaskGeometry
vector < vector <unsigned char> > get_detector_mask(int quadID, int detNo, TanMaskGeometry tanMask, int &status){
    int row, col;
    int startx, starty;
    int ipixx=0, ipixy=0; //counter variable
    vector < vector <unsigned char> > detMask;
    vector <unsigned char> tempMask;
    
    //Initializing vector
    for(ipixx=0; ipixx<DET_MASKSIZE; ipixx++){
        tempMask.push_back(0);
    }
    for(ipixy=0; ipixy<DET_MASKSIZE; ipixy++){
        detMask.push_back(tempMask);
    }
    //Vectors initialized
    
    if(quadID<0 || quadID>3){
        LOG(ERROR)<<"[ERROR] Quadrant ID should be either 0 or 1 or 2 or 3";
        status = EXIT_FAILURE;
        exit(EXIT_FAILURE);
    }
    if(detNo<0 || detNo>15){
        LOG(ERROR)<<"[ERROR] Wrong detector ID value. Detector ID should lie between 0 to 15";
        status = EXIT_FAILURE;
        exit(EXIT_FAILURE);
    }
    col = detNo%NUM_DET_PER_QUAD_PER_COL;
    row = detNo/NUM_DET_PER_QUAD_PER_ROW;
    
    startx = col*(DET_MASKSIZE + NPIXELS_GAP);
    starty = row*(DET_MASKSIZE + NPIXELS_GAP);
    
    DLOG(INFO) << "Extracting mask for detector " << detNo << " of quadrant " << quadID;
    for(ipixy=starty; ipixy<(starty+DET_MASKSIZE); ipixy++){
        for(ipixx=startx; ipixx<(startx+DET_MASKSIZE); ipixx++){
            if (quadID==0){
                detMask[ipixy-starty][ipixx-startx] = tanMask.q0UnMask[ipixy][ipixx];
            }
            else if (quadID==1){
                detMask[ipixy-starty][ipixx-startx] = tanMask.q1UnMask[ipixy][ipixx];
            }
            else if (quadID==2){
                detMask[ipixy-starty][ipixx-startx] = tanMask.q2UnMask[ipixy][ipixx];
            }
            else if (quadID==3){
                detMask[ipixy-starty][ipixx-startx] = tanMask.q3UnMask[ipixy][ipixx];
            }
        }
    }
        
    return detMask;
}

vector < vector <unsigned char> > get_shifted_detector_mask(int quadid, int detNo, TanMaskGeometry tanMask, float dX, float dY, int& status){
    vector <vector <unsigned char> > detMask;
    vector <vector <unsigned char> > shiftedDetMask;
    vector <unsigned char> tempMask;
    int rowstart, rowend, colstart, colend; //for shifted mask
    int xstart, ystart;
    int i=0, j=0, x=0 ,y=0; //
    int shiftPixX=0;
    int shiftPixY=0;
    //Getting detector mask for a particular quadrant;
    detMask = get_detector_mask(quadid, detNo, tanMask, status);
    if(status){
        LOG(ERROR) << "Error in getting detector mask for detector " << detNo << " of quadrant " << quadid;
        exit(EXIT_FAILURE);
    }
    if(dX==0 && dY==0){
        return detMask;
    }
    
    shiftPixX = (int) abs(dX/PIX_SIZ_MASK);
    shiftPixY = (int) abs(dY/PIX_SIZ_MASK);
    
    DLOG(INFO) << "Shift in pixels: shiftX=" << shiftPixX << " shiftY=" << shiftPixY;
    
    //Storing detector mask in shiftedDetMask without applying shifts
    for(i=0; i<DET_MASKSIZE;i++){
        tempMask.push_back(0);
    }
    for(i=0; i<DET_MASKSIZE; i++){
        shiftedDetMask.push_back(tempMask);
    }


    if (dX > 0) {
        colstart = 0;
        colend = DET_MASKSIZE - shiftPixX;
        xstart = shiftPixX;
    } else {
        colstart = shiftPixX;
        colend = DET_MASKSIZE;
        xstart = 0;
    }

    if (dY > 0) {
        rowstart = 0;
        rowend = DET_MASKSIZE - shiftPixY;
        ystart = shiftPixY;
    } else {
        rowstart = shiftPixY;
        rowend = DET_MASKSIZE;
        ystart = 0;
    }    
    
  for(i=rowstart,y=ystart;i<rowend;i++,y++){
        for(j=colstart,x=xstart;j<colend;j++,x++){
            shiftedDetMask[i][j]=detMask[y][x];
            if(i>DET_MASKSIZE || j>DET_MASKSIZE){
                LOG(ERROR)<<"Index out of range for shifted mask";
                status = EXIT_FAILURE;
               exit(EXIT_FAILURE);
            }
             if(y>DET_MASKSIZE || x>DET_MASKSIZE){
                LOG(ERROR)<<"Index out of range for  mask";
                status = EXIT_FAILURE;
                exit(EXIT_FAILURE);
            }
            
        }
    }

    return shiftedDetMask;
}



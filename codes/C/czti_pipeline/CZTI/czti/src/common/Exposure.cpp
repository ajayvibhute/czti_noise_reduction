#include "Exposure.h"

using namespace std;

//SHADOW CLASS

Shadow::Shadow(string cameraGeomfile,float sourceFlux, float thetaXr, float thetaYr,
        int nBins, float energyStart, float energyEnd, int includeMask, int maskOversampling) {
    this->cameraGeomfile = cameraGeomfile;
    this->sourceFlux = sourceFlux;
    this->thetaXr=thetaXr;
    this->thetaYr=thetaYr;
    this->nBins = nBins;
    this->energyStart = energyStart;
    this->energyEnd = energyEnd;
    this->includeMask = includeMask;
    this->maskOversampling = maskOversampling;
    
}

int Shadow::display() {
    LOG(INFO) << "---------------------------------------------------------------------";
    LOG(INFO) << "                       SHADOW COMPUTATION PARAMETERS                 ";
    LOG(INFO) << "---------------------------------------------------------------------";
    LOG(INFO) << "Camera Geometry File                        : " << cameraGeomfile;
    LOG(INFO) << "Thetax, Thetay (in radians)                 : " << thetaXr << ", " << thetaYr;
    LOG(INFO) << "Thetax, Thetay (in degrees)                 : " << thetaXr*TODEG << ", " << thetaYr*TODEG;
    LOG(INFO) << "Source Flux                                 : " << sourceFlux;
    LOG(INFO) << "Energy Start                                : " << energyStart;
    LOG(INFO) << "Energy End                                  : " << energyEnd;
    LOG(INFO) << "Number of energy bins                       : " << nBins;
    LOG(INFO) << "Whether to include mask?                    : " << includeMask;
    LOG(INFO) << "Oversampling factor (mask transmission)     : " << maskOversampling;
    LOG(INFO) << "---------------------------------------------------------------------";
    
}

int Shadow::check_shadow_par() {
    int flag=0;
    if(cameraGeomfile==""){
        LOG(ERROR) << "Enter valid camera geometry file name.";
        flag=1;
        return flag;
    } else if( sourceFlux <= 0.0){
        LOG(ERROR) << "Source flux should be greater than 0.";
        flag=1;
        return flag;
    } else if(energyStart<10.0 || energyEnd>100.0 || energyStart>energyEnd){
        LOG(ERROR) << "Valid energy range is 10.0 to 100.0 keV.";
        flag=1;
        return flag;
    } else if(nBins<1 || nBins>20){
        LOG(ERROR) << "Number of energy bins should be between 1 & 20";
        flag=1; return flag;
    } else if(includeMask==1){
        if(maskOversampling<2 || maskOversampling>20){
            LOG(ERROR) << "Mask oversampling should be in the range 2 to 20.";
            flag=1; return flag;
        }
    }
    return flag;
}

int Shadow::compute_shadow(vector <vector <float> > &fullnormalizedEffarea,
        vector <vector <unsigned char> > &fullfinermask){
    int status=0;
    TanMaskGeometry msk; //to get full uncompressed mask
    ExposureTable expTable; //to store openfraction values 
    float thetaXd, thetaYd; //thetax and thetay in degrees
    float thetar;
    int idetx=0, idety=0, ipixx=0, ipixy=0, ix=0, iy=0;
    vector <float> tempShadow;
    float shadowPix=0.0;
    
    //Checking parameter values
    if(check_shadow_par()){
        LOG(ERROR) << "Error in shadow parameters supplied";
        return EXIT_FAILURE;
    }
    //initializing shadow vector
    tempShadow.resize(NPIXELS_X, 0.0);
    shadowTab.shadow2D.resize(NPIXELS_Y, tempShadow);

    //converting thetax and thetay into degrees
    thetaXd = thetaXr * TODEG;
    thetaYd = thetaYr * TODEG;
    
    //calculating theta
    thetar = atan(tan(thetaXr) * tan(thetaXr) + tan(thetaYr) * tan(thetaYr));
    
    if (nBins > 1) {
        if (calculate_average_pixexposure(thetaXd, thetaYd, cameraGeomfile, fullfinermask,
                expTable, includeMask, maskOversampling, nBins, energyStart, energyEnd)) {
            LOG(ERROR) << "Error in evaluating exposure map.";
            return EXIT_FAILURE;
        }
    } else if (nBins == 1) {
        if (calculate_average_pixexposure(thetaXd, thetaYd, cameraGeomfile, fullfinermask,
                expTable, includeMask, maskOversampling, nBins, energyStart)) {
            LOG(ERROR) << "Error in evaluating exposure map.";
            return EXIT_FAILURE;
        }
    }

    if(expTable.generate_full_exposure_image()){
        LOG(ERROR) << "Error in generating full exposure image.";
        return EXIT_FAILURE;
    }
    
    for (idety = 0; idety < NDET_PER_COL; idety++) {
        for (idetx = 0; idetx < NDET_PER_ROW; idetx++) {
            for (ipixy = 0; ipixy < PIXELS_PER_COL_DET; ipixy++) {
                for (ipixx = 0; ipixx < PIXELS_PER_ROW_DET; ipixx++) {
                    ix =idetx*PIXELS_PER_ROW_DET + ipixx;
                    iy =idety*PIXELS_PER_COL_DET + ipixy;
                    
                    //Calculating shadow value for each pixel
                    shadowPix = sourceFlux * expTable.fullExposure[iy][ix] * fullnormalizedEffarea[iy][ix];
                    shadowTab.shadow2D[iy][ix] = shadowPix;
                    shadowTab.set_shadow(idetx, idety, ipixx, ipixy, ix, iy, shadowTab.shadow2D[iy][ix]);
                }
            }

        }
    }
    //setting information parameters of shadow table
    shadowTab.thetaXr = this->thetaXr;
    shadowTab.thetaYr = this->thetaYr;
    shadowTab.energyMin = this->energyStart;
    shadowTab.energyMax = this->energyEnd;
    shadowTab.nBins = this->nBins;
    shadowTab.sourceFlux = this->sourceFlux;
    

    return EXIT_SUCCESS;
}
//SHADOW CLASS END

//INDEPENDENT FUNCTIONS
int calculate_pixexposure(float ekev, float thetax, float thetay, 
        string cameraGeomFileName, vector < vector <float> > &ofmaskFull, int includeMask,
        ExposureTable &expTable){
    float offset; //offset for thetax to avoid degeneracy of surface identity
    float tt; //temporary variable
    int isurf, nsurf; //surface counter variable, number of camera surfaces
    int idetx, idety, ipixx, ipixy, ix, iy, imaskx, imasky; //counter variables
    long nrows; //number of rows in fits table
    double tradx, trady; //thetax and thetay in radians
    double ctx, cty, stx, sty, tgx, tgy; //cos(thetax), cos(thetay), sin, tan.... 
    double dfac; //secant of the net angle from the CZTI axis;
    double vl, vm, vn; //direction cosines
    double alphaAl, alphaTa, alphaCZT, alphaPCB; //absorption coefficients of Al, Ta, CZT & PCB
    double tauAl, tauTa, tauCZT, tauPCB; //optical depth of corresponding materials
    double cztabs, trTam; //czt absorption fraction, Transmission fraction of tantalummask
    double xmoff, ymoff; //offset w.r.t. the detector plane
    double pr; //temporary variable
    vector<double> proj; // cosine of angle between the vector and the normal to a surface
    vector<double> strans; // strans is the transmission through the surface at the given energy and angle
    int nx, ny; // to subdivide into nx x ny subpixles
    float x0, y0; //bottom-left coordinate of a detector
    float xpix0, ypix0; //bottom-left coordinate of a pixel
    float xpixw, ypixw; //length & breadth of a pixel
    float dx, dy; //edge length of each sub-pixel
    double x,y; //bottom-left coordinate of each sub-pixel;
    double Ageopix; //Geometric area of a pixel
    double trans;  //net transmission to an area bin through all intervening surfaces.
    double farea;
    double eta; //length of the vector from teh current area element to the surface.
    double xint, yint, zint; //co-ordinates of point of intersection
    double xmin, xmax, ymin, ymax; 
    float bin; //to store area of each sub-pixel
    double ofmask;
    double pixarea, pixexposure;
    double weight; //weights without taking into account renormalization offset values.
    //to store information in camera geometry CALDB file
    CztiCameraGeometry camGeometry; 
    vector <unsigned short> indexSurf;
    vector <float> zmin, zmax;
    vector <float> xminz, xminc, xmaxz, xmaxc;
    vector <float> yminz, yminc, ymaxz, ymaxc;
    vector <float> xnp, ynp, znp, dnp;
    vector <float> thickAl, thickTa, openfrac;
    
    //to store information regarding czti detector geometry
    CztDetectorGeometry cztDetGeometry;
    vector <float> xdet0; //x coordinate of 8 detector locations (cm)
    vector <float> ydet0; //y coordinate of 8 detector locations (cm)
    vector <float> pixx0; //x coordinate of 16 pixel locations (cm)
    vector <float> pixy0; //y coordinate of 16 pixel locations (cm)
    vector <float> pixxw; //width of pixel in x direction
    vector <float> pixyw; //width of pixel in y direction
    
    //absorption coefficients of materials used
    alphaAl=absco_Al(ekev);
    alphaTa=absco_Ta(ekev);
    alphaCZT=absco_CZT(ekev);
    alphaPCB=absco_pcb(ekev);
    
    //limit angle range to avoid zero divide;
    if(abs(90.0-abs(thetax)) < 1.0){
        thetax=89.50;
    }
    if(abs(90.0-abs(thetay)) < 1.0){
        thetay=89.50;
    }
    
    //Avoid degeneracy of surface identity at azimuth diagonals by applying a small offset to thetax
    offset=0.01;
    if (thetax == thetay) {
        if (thetax != 0.0) {
            tt = abs(thetax) - abs(offset);
            if (thetax < 0.0) {
                thetax = -abs(tt);
            } else {
                thetax = abs(tt);
            }

        }
    }
    
    //calculating trigonometric ratios
    tradx = thetax * M_PI / 180;
    trady = thetay * M_PI / 180;
    ctx = cos(tradx);
    stx = sin(tradx);
    tgx = stx / ctx;
    cty = cos(trady);
    sty = sin(trady);
    tgy = sty / cty;
    
    //dfac is the secant of the net angle from the CZTI axis
    dfac = sqrt(1.0+tgy*tgy+tgx*tgx);
    
    //direction cosines vl, vm, vn of the vector pointing in the direction 
    //thetax, thetay from normal incidence.
    vl=tgx/dfac;
    vm=tgy/dfac;
    vn=1.0/dfac;
    
    // effective CZT detector thickness
    tauCZT=abs(CZTDET_THICKNESS*alphaCZT*dfac);
    cztabs=1.0-exp(-tauCZT);
    
    // PCB layer of thickness 0.015 cm on CZT detector dimnishes the radiation 
    // reaching the detector, and hence its effective absorption
    tauPCB=abs(PCB_THICKNESS*alphaPCB*dfac);
    cztabs=cztabs*exp(-tauPCB);
    
    cztabs=1.0;
    // Mask plate is made of Tantalum of thickness 0.05 cm. trTam is the 
    // transmission fraction through a closed portion of the mask plate.
    trTam=exp(-1*TANMASK_THICKNESS*alphaTa*dfac);
    
    //offset w.r.t. the detector plane of the ray intersection point with mask plate.
    xmoff=tgx*ZMASK;
    ymoff=tgy*ZMASK;
    
    //Reading and storing information from camera geometry file
    //cameraGeomFileName = caldb_full_path_generator(cameraGeomFileName);
    if(camGeometry.read_camera_geometry_file(cameraGeomFileName)){
        LOG(ERROR) << "Error in reading camera geometry CALDB file: "<< cameraGeomFileName;
        return(EXIT_FAILURE);
    }
    indexSurf = camGeometry.get_indexSurf();
    zmin = camGeometry.get_zmin();
    zmax = camGeometry.get_zmax();
    xmaxz = camGeometry.get_xmaxz();
    xmaxc = camGeometry.get_xmaxc();
    ymaxz = camGeometry.get_ymaxz();
    ymaxc = camGeometry.get_ymaxc();
    xminz = camGeometry.get_xminz();
    xminc = camGeometry.get_xminc();
    yminz = camGeometry.get_yminz();
    yminc = camGeometry.get_yminc();
    xnp   = camGeometry.get_xnp();
    ynp   = camGeometry.get_ynp();
    znp   = camGeometry.get_znp();
    dnp   = camGeometry.get_dnp();
    thickAl = camGeometry.get_thickAl();
    thickTa = camGeometry.get_thickTa();
    openfrac = camGeometry.get_openfrac();
    
    //Number of camera surfaces
    nsurf = indexSurf.size();
    
    
    //initializing proj
    proj.clear();
    strans.clear();
    for(isurf=0; isurf < nsurf; isurf++){
        proj.push_back(0.0);
        strans.push_back(0.0);
    }
    for (isurf = 0; isurf < nsurf; isurf++) {
        
        //Pre compute and store transmission through each surface
        
        //proj is the cosine of the angle between the vector and the normal to a surface
        pr=xnp[isurf]*vl + ynp[isurf]*vm + znp[isurf]*vn;
        //avoid zero divide
        if(pr==0.0) { pr=1e-9;}
        proj[isurf]=pr;
        
        //compute optical depth of the Aluminium and Tantalum components
        tauAl=abs(alphaAl*thickAl[isurf]/pr);
        tauTa=abs(alphaTa*thickTa[isurf]/pr);
        
        //from optical depth compute the transmission fraction. Note that some surfaces
        //may have material with partial cover (e.g. coded mask), which is parameterized
        //by openfrac
        
        strans[isurf] = exp(-tauAl)*(openfrac[isurf] + (1.0-openfrac[isurf])*exp(-tauTa));    
    }
    
    //Initializing detector geometry
    if(cztDetGeometry.init_detectors()){
        LOG(ERROR) << "Error in initializing CZT detectors";
    }
    xdet0 = cztDetGeometry.get_xdet0();
    ydet0 = cztDetGeometry.get_ydet0();
    
    //each pixel is subdivided into nx x ny sub pixels for exposure computation
    nx=4;
    ny=4;
    
    //loop over detectors: 8x8 array
    for(idetx=0; idetx<8; idetx++){
        x0=xdet0[idetx];
        for(idety=0; idety<8; idety++){
            y0=ydet0[idety];
            
            cztDetGeometry.define_pixel_coordinates(x0,y0);
            pixx0 = cztDetGeometry.get_pixx0();
            pixy0 = cztDetGeometry.get_pixy0();
            pixxw = cztDetGeometry.get_pixxw();
            pixyw = cztDetGeometry.get_pixyw();
            
            for(ipixx=0; ipixx<16; ipixx++){
                for(ipixy=0; ipixy<16; ipixy++){
                    //bottom-left coordinates of each pixel
                    xpix0=pixx0[ipixx];
                    ypix0=pixy0[ipixy];
                    //edge length and breadth of each pixel
                    xpixw = pixxw[ipixx];
                    ypixw = pixyw[ipixy];
                    
                    //edge length and breadth of each sub-pixel
                    dx=pixxw[ipixx]/(1.0*nx);
                    dy=pixyw[ipixy]/(1.0*ny);
                    
                    //Geometric area of a pixel
                    Ageopix=xpixw*ypixw; //Geometric area of a pixel
                    
                    //CORE ALGORITHM OF EFFECTIVE AREA COMPUTATION
                    
                    //on each pixel, integrate over detector surface area, adding effective
                    //contribution from one area bin at a time
                    farea=0.0;
                    //Geometric area of each sub-pixel
                    bin=dx*dy;
                    //bottom-left coordinates of each sub-pixel
                    for(ix=0; ix<nx; ix++){
                        x=xpix0+(1.0*ix + 0.5)*dx;
                        for(iy=0; iy<ny; iy++){
                            y=ypix0+(1.0*iy + 0.5)*dy;
                            
                            //Here, (x,y) is the center of the current area bin
                            
                            //trans is the net transmission to an area bin through all the 
                            //intervening surfaces. It is a product of transmission fractions
                            //of all the intervening surfaces.
                            
                            trans=1.0;
                            
                            //loop over the surfaces
                            for(isurf=0; isurf<nsurf; isurf++){
                                //eta is the length of the vector from the current area 
                                //element to the surface.
                                eta=-(xnp[isurf]*x+ynp[isurf]*y-dnp[isurf])/proj[isurf];
                                
                                //coordinates of the point of intersection between the vector
                                //and the surface
                                xint=x+vl*eta; //X intersection
                                yint=y+vl*eta; //Y intersection
                                zint=vn*eta; //Z intersection
                                
                                //Check if intersection falls between the surface boundaries zmin,
                                //zmax, ymin, ymax, xmin, xmax. For slanted surfaces the x-y limits
                                //depend on z, and are computed from the equations below.
                                
                                if((zint>=zmin[isurf]) && (zint<=zmax[isurf])) {
                                    xmin = xminz[isurf] * zint + xminc[isurf];
                                    xmax = xmaxz[isurf] * zint + xmaxc[isurf];
                                    ymin = yminz[isurf] * zint + yminc[isurf];
                                    ymax = ymaxz[isurf] * zint + ymaxc[isurf];
                                    
                                    //if surface is intersected within its boundary then take into 
                                    //account the fractional transmission through it.
                                    if(xint>=xmin && xint<=xmax && yint>=ymin && yint<=ymax){
                                        trans=trans*strans[isurf];
                                    }
                                }
                                
                            }//END LOOP OVER SURFACES
                            
                            //Mask transmission is ignored (set to unity)
                            /*--------------------------------------------
                             //Compute transmission through the Coded mask plate
                             xint=x+xmoff;
                             yint=y+ymoff;
                             trans=trans*transmask(xint, yint, trTam);
                            --------------------------------------------- */
                            
                            //increment the effective area by adding transmission 
                            //fraction times the geometric are of the bin
                            farea=bin*trans+farea;
                            //DLOG(INFO) << "farea for (ipixx, ipixy, ix, iy)" << ipixx << " "<<
                            //        ipixy << " " << ix << " " << iy << " :" << farea; 
                        }//END LOOP OVER Y-POSITION OF AREA BIN ON A PIXEL   
                    }//END LOOP OVER X-POSITION OF AREA BIN ON A PIXEL
                    
                    /*
                    //----------------------------------------------------------------
                    //multiply farea by cosine of the angle from the CZTI axis to take
                    //care of the geometric projection (1/dfac).
                    //Also multiply by the absorption efficiency of the CZT detector
                    //at this angle & energy
                    
                            pixarea=farea*cztabs/dfac;
                    //----------------------------------------------------------------
                     */
                    //Masktransmission
                    if (includeMask == 1) {
                        imaskx = idetx * NO_PIX_X_PER_DET + ipixx;
                        imasky = idety * NO_PIX_Y_PER_DET + ipixy;
                        ofmask = (double) ofmaskFull[imasky][imaskx];;
                        farea = farea * ((ofmask) + (1-ofmask)*trTam);
                    }
                    pixarea=farea*cztabs;
                    //END OF CORE EFFECTIVE AREA ALGORITHM
                    
                    //CORRECTION OF RELATIVE QE & THRESHOLD REMOVED [FORTRAN]
                    //<MODIFY IT TO USE IT>
                    /*---------------------------------------------------------
                        pixarea=relqe(idetx,idety,ipixx,ipixy)*pixarea
                        if (thresh(idetx,idety,ipixx,ipixy).gt.ekev) then
                            pixarea=0.0
                        endif

                     ---------------------------------------------------------*/
                    pixexposure=pixarea/Ageopix;
                    weight = (2.0*pixexposure)-1.0;
                    //LOG(INFO) << imaskx << " " << imasky << " " << setprecision(10) << pixexposure << " " << weight;
                    
                    expTable.set_exposure((unsigned char)idetx, (unsigned char)idety, 
                            (unsigned char)ipixx, (unsigned char)ipixy, 
                            (float) pixexposure, (float) weight);
                }//END LOOP OVER Y PIXEL NUMBER
            }//END LOOP OVER X PIXEL NUMBER
        }//END LOOP OVER DETECTOR ARRAY COLUMN (y)
    }//END LOOP OVER DETECTORY ARRAY ROW (x)
    

    return EXIT_SUCCESS;
}

int create_exposure_array(float thetaxd, float thetayd, string cameraGeomFilename, 
        vector <vector <unsigned char> > &fullfinermask, ExposureTable &expTable, int includeMask, int maskOversampling, 
        string eboundsFilename, float minEnergy, float maxEnergy, int nchannels){
    int status = 0;
    int ichannel= 0; //index based on channel
    int i=0;
    float tempN = 0.0;
    string energyRange="";
    float minE=0.0; //minimum energy
    float maxE=0.0; //maximum energy
    TransMask transmask;
    Ebounds ebounds; //to read and store ebounds data
    bool efile=FALSE;
    long numchannels=0; //number of channels
    CztDetectorGeometry cztdet; //to store coordinates of detector pixels
    TanMaskGeometry msk; //to get full uncompressed mask
    vector <vector <float> > fullofmask; //full mask open fraction with transmission
    ExposureTable tempexposure; //temporary vector to find average exposure of all energy bins.
    vector <vector <float> > openfractionArray;
    vector <vector <float> > weightsArray; //weights without renormalization offset
    vector <float> energiesArray; //energy values (in keV) for which openfractionArray is calculated
    vector <unsigned short> PIsArray; //PI bins for which openfractionArray is generated.
    vector <float> tempof;
    float binsize =0.0;
    float ekev = 0.0;
    unsigned int channel=0;
    
    if(eboundsFilename=="") {
        //Energy validations
        if (minEnergy < 1.0) {
            LOG(ERROR) << "Minimum energy cannot be less than 1.0 keV.";
            return EXIT_FAILURE;
        }
        if (minEnergy > maxEnergy) {
            LOG(ERROR) << "Minimum energy cannot be less than maximum energy.";
            return EXIT_FAILURE;
        }
        if(nchannels<0) {
            LOG(ERROR) << "Number of channels cannot be negative";
            return EXIT_FAILURE;
        }
        numchannels = nchannels;
        energyRange = itoa(minEnergy) + "-" + itoa(maxEnergy);
        efile = FALSE;
    }
    else if(FileExists((char*) eboundsFilename.c_str())){
        if(ebounds.read_ebounds_file(eboundsFilename)){
            LOG(ERROR) << "Error in reading Ebounds file: " << eboundsFilename;
            return EXIT_FAILURE;
        }
        numchannels = ebounds.get_nrows(); 
        minE = ebounds.get_energy(0);
        maxE = ebounds.get_energy(numchannels-1);
        energyRange = itoa(minE)+"-"+itoa(maxE);
        efile = TRUE;
    } else {
        LOG(ERROR) << "Ebounds file specified by user does not exist: " << eboundsFilename;
        return EXIT_FAILURE;
    }
    
    tempof.resize(numchannels,0.0);
    openfractionArray.resize(XSIZE*YSIZE, tempof);
    weightsArray.resize(XSIZE*YSIZE, tempof);
    energiesArray.resize(numchannels, 0.0);
    PIsArray.resize(numchannels, 0);
    
    //Getting detector geometry
    if (cztdet.init_detectors()) {
        LOG(INFO) << "Error in getting detector coordinates (physical value).";
        return EXIT_FAILURE;
    }

    //setting thetax and thetay for mask fraction calculation
    transmask.set_thetaXY((double) thetaxd, (double) thetayd);

    //storing fullfinermask as a variable of transmask object
    transmask.set_fullUnMask(fullfinermask);
    //Calculating mask open fraction

    //Calculating mask open fraction
    if (transmask.compute_mask_openfraction(cztdet, maskOversampling)) {
        LOG(ERROR) << "Error in calculating mask open fraction.";
        return EXIT_FAILURE;
    }

    if (transmask.generate_full_openfraction_image()) {
        LOG(ERROR) << "Error in generating mask openfraction image.";
        return EXIT_FAILURE;
    }
    fullofmask.clear();
    fullofmask = transmask.get_full_openfraction();

    for (ichannel=0; ichannel<numchannels; ichannel++) {

        tempexposure.reset();
        if(efile==TRUE){
			//LOG(INFO)<<"Getting from ebounds";
            try{
                ekev = ebounds.get_energy(ichannel);
                channel = ebounds.get_channel(ichannel);
            } catch (string s){
                LOG(ERROR) << s;
                return EXIT_FAILURE;
            }
//			LOG(INFO)<<"Getting from ebounds row "<<ichannel<<" channel "<<channel<< " Energy "<<ekev;
            
        } else if(efile==FALSE){
            binsize = (maxEnergy-minEnergy)/((float)numchannels);
            ekev = minEnergy + (ichannel+0.5)* binsize;
            channel=ichannel;
        }
        
        DLOG(INFO) << "Evaluating open fraction for energy :" << ekev << "keV";

        //Calculating pixel exposure i.e. openfraction for each pixel 
        if (calculate_pixexposure(ekev, thetaxd, thetayd, cameraGeomFilename, fullofmask,
                includeMask, tempexposure)) {
            LOG(ERROR) << "Error in evaluating exposure map.";
            return EXIT_FAILURE;
        }
        
        energiesArray[(int) ichannel] = ekev;
        PIsArray[(int) ichannel]= channel;
        for(i=0; i<XSIZE*YSIZE; i++){
            openfractionArray[i][(int)ichannel] = tempexposure.openfrac[i];
            weightsArray[i][(int)ichannel] = tempexposure.weights[i];
            
        }
    }

    //Copying values to exptable;
    expTable.reset();
    expTable = tempexposure;
    expTable.set_openFracArray(openfractionArray);
    expTable.set_weightsArray(weightsArray);
    expTable.PIs = PIsArray;
    expTable.thetaxd = thetaxd;
    expTable.thetayd = thetayd;
    expTable.energyRange = energyRange; 
    expTable.nchannels = numchannels;
    expTable.energies = energiesArray;
    expTable.eboundFlag = efile?true:false;
    //clearing temporary vectors
    tempexposure.reset();
    
    return status;
}

int calculate_average_pixexposure(float thetax, float thetay, string cameraGeomFilename, 
        vector <vector <unsigned char> > &fullfinermask, ExposureTable &expTable, int includeMask, int maskOversampling, int nbins, 
        float minEnergy, float maxEnergy){
    int status=0;
    int iEnergy=0; //index based on energy
    float tempN=0.0;
    TransMask transmask;
    CztDetectorGeometry cztdet; //to store coordinates of detector pixels
    TanMaskGeometry msk; //to get full uncompressed mask
    vector <vector <float> > fullofmask; //full mask open fraction with transmission
    ExposureTable tempexposure; //temporary vector to find average exposure of all energy bins.
    vector <float> tempopenfraction;
    float ekev=0.0;
    string energyRange;

    //Getting detector geometry
    if (cztdet.init_detectors()) {
        LOG(INFO) << "Error in getting detector coordinates (physical value).";
        return EXIT_FAILURE;
    }
    
    //setting thetax and thetay for mask fraction calculation
    transmask.set_thetaXY((double) thetax, (double) thetay);
    
    //storing fullfinermask as a variable of transmask object
    transmask.set_fullUnMask(fullfinermask);
    //Calculating mask open fraction

    //Resizing tempopenfraction and initializing it with value 0.0
    tempopenfraction.resize(NO_PIX_ALL_QUADS, 0.0);
    
    //Calculating mask open fraction
    if (transmask.compute_mask_openfraction(cztdet, maskOversampling)) {
        LOG(ERROR) << "Error in calculating mask open fraction.";
        return EXIT_FAILURE;
    }

    if (transmask.generate_full_openfraction_image()) {
        LOG(ERROR) << "Error in generating mask openfraction image.";
        return EXIT_FAILURE;
    }
    fullofmask.clear();
    fullofmask = transmask.get_full_openfraction(); 
    
    for (iEnergy = 0; iEnergy < nbins; iEnergy++) {

        tempexposure.reset();
        if(nbins>1) {
        ekev=minEnergy + ((maxEnergy-minEnergy)/(nbins-1))*iEnergy;
        }
        else if(nbins==1){
            ekev = minEnergy;
        }
        LOG(INFO) << "Evaluating openfraction for energy :" << ekev << "keV";
        
        //Calculating pixel exposure i.e. openfraction for each pixel 
        if (calculate_pixexposure(ekev, thetax, thetay, cameraGeomFilename, fullofmask, 
                includeMask,tempexposure)) {
            LOG(ERROR) << "Error in evaluating exposure map.";
            return EXIT_FAILURE;
        }
        
        //Adding openfraction (of each pixel) for all energy bins. 
        transform(tempopenfraction.begin(), tempopenfraction.end(),
                tempexposure.openfrac.begin(), tempopenfraction.begin(), plus<float>());
    }

    
    tempN=1.0/nbins;
    //Taking average of all openfraction values
    LOG(INFO) << "Taking average of openfraction values for number of energy bins " << nbins;
    transform(tempopenfraction.begin(), tempopenfraction.end(),
            tempopenfraction.begin(), bind1st(multiplies<float>(),tempN));
    
    if(nbins>1)
        energyRange = itoa(minEnergy) + "-" + itoa(maxEnergy);
    else if(nbins==1)
        energyRange = itoa(ekev);
    //Copying values to exptable;
    expTable.reset();
    expTable = tempexposure;
    expTable.set_openfrac(tempopenfraction); //average openfraction value;
    expTable.thetaxd = thetax;
    expTable.thetayd = thetay;
    expTable.energyRange = energyRange; 

    
    //clearing temporary vectors
    tempexposure.reset();
    
    return EXIT_SUCCESS;
}

int calculate_renormalized_weights(ExposureTable &expTable, string badpixFilename, int badpixThreshold, 
        vector <int> quadsToProcess, string effareaFilename, bool includeCztPcbAbsorption,string evtfilename){
    int status=0;
    int i=0;
    int iquad=0, ipix=0, ienergy=0;
    int pixx=0; //location of pixel in a quadrant (0-63))
    int pixy=0;
    int locx=0; // location of pixel in all quadrants (0-127)
    int locy=0;
    Badpix badpix; 
    int pixquad=0; //quadrant in which pixel is present
    vector < vector <unsigned char> > badpixMap; //badpixel map for all 128x128 pixels
    vector <float> zeroWeights;
    vector <EffArea> effarea; //to read effective area file
    vector <float> renormOffsets; //renormalization offsets at all energies specified by user
    vector <float> areas; //sum of area of all pixels at particular energies
    char extname[FLEN_VALUE];
    vector <float> effectiveArea; //effective area of a particular quadrant at all energies for which
                                  //weights have to evaluated.
    vector <float> tempEnergies;
    effarea.resize(4); //to store effecive area values for all 4 quadrants
    zeroWeights.resize(expTable.nchannels, 0.0);

    double exposurefrac_locxy[128][128];
    double exposure_fraction[4][4096];
    int detid,pixid;


    if(badpix.read_badpix_file(badpixFilename)){
        LOG(ERROR) << "Error in reading bad pixel file " << badpixFilename;
        return EXIT_FAILURE;
    }
    
    badpixMap = badpix.get_badpix_map();
    if(badpixMap.empty()){
        LOG(INFO) << "Creating badpix map.";
        try {badpix.create_badPixMap();}
        catch(string s){
            LOG(ERROR)<< s;
            return EXIT_FAILURE;
        }
        badpixMap = badpix.get_badpix_map();
    }
    
    expTable.badpixFlag.resize(NO_PIX_ALL_QUADS, 0);
    //setting weights of exposure table to 0 for all badpixels as defined by badpixThreshold
    //This is done irrespective of quadrants to be processed
    for (i = 0; i < NUMQUAD; i++) {
        iquad = i;
        for (pixx = 0; pixx < XPIX_QUAD; pixx++) {
            for (pixy = 0; pixy < YPIX_QUAD; pixy++) {
                generate_locx_locy(pixx, iquad, pixy, locx, locy);
                for (ipix = 0; ipix < expTable.get_nrows(); ipix++) {
                    if (expTable.locX[ipix] == locx && expTable.locY[ipix] == locy) {
                        expTable.badpixFlag[ipix]= badpixMap[locy][locx];
                        if ((int) badpixMap[locy][locx] > badpixThreshold) {
                            expTable.weightsArray[ipix] = zeroWeights;
                        } //END IF 
                    }//END IF ON BADPIXEL THRESHOLD CONDITION
                }//END LOOP ON IPIX (128X128) 
            }//END LOOP ON PIXY (SINGLE QUADRANT)
        }//END LOOP ON PIXX (SINGLE QUADRANT)
    }//END LOOP ON IQUAD (QUADRANTS TO BE PROCESSED)
    
    //Reading effective area file
    for (i = 0; i < quadsToProcess.size(); i++) {
        iquad = quadsToProcess[i];
        quadToHDU(iquad, extname);
        if ((effarea[iquad]).read_effarea_file(effareaFilename, (string) extname)) {
            LOG(ERROR) << "Error in reading effective area file: " << effareaFilename <<
                    " for quadrant " << iquad;
            return EXIT_FAILURE;
        }
    }
    //Effective area file for all quadrants read.

    //Read the exposure fractions of each pixel

    fitsfile *fevt;

    fits_open_file(&fevt, (char *)evtfilename.c_str(), READONLY, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in opening event file";
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fevt,BINARY_TBL,"EXPOSURE",0,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to EXPOSURE extension of fits file";
        return (EXIT_FAILURE);
    }
    
    long nrows;
    int detx,dety;

    fits_get_num_rows(fevt, &nrows, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in getting numrows";
        return (EXIT_FAILURE);
    }

    for(i = 0; i < quadsToProcess.size(); i++)
    {

        char colname[20];
        int exp_col;
        sprintf(colname,"EXPOSURE_Q%d",quadsToProcess[i]);

        fits_get_colnum(fevt,CASEINSEN,colname,&exp_col,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Not finding exposure fraction in the event file";
            return (EXIT_FAILURE);
        }

        fits_read_col(fevt, TDOUBLE, exp_col, 1, 1, nrows, NULL, exposure_fraction[i],NULL, &status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in moving to EXPOSURE extension of fits file";
            return (EXIT_FAILURE);
        }
    

        //Convert the exposure into locx locy format
        for(detid=0;detid<16;detid++)
        {
            for(pixid=0;pixid<256;pixid++)
            {
                generate_locx_locy(detid, pixid,quadsToProcess[i], detx, dety, locx, locy);
                exposurefrac_locxy[locx][locy]=exposure_fraction[quadsToProcess[i]][detid*256+pixid];   
            }    
        }
    }

    fits_close_file(fevt,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }
 
    //Done reading exposure fractions

       
    //Calculating renormalization offsets 
    renormOffsets.resize(expTable.energies.size(), 0.0);
    areas.resize(expTable.energies.size(), 0.0);
    tempEnergies.resize(expTable.energies.size(), 30.0);
   
    for (i = 0; i < expTable.get_nrows(); i++){
        if(is_pix_present_in_quads(expTable.locX[i], expTable.locY[i], quadsToProcess,
                pixquad, pixx, pixy)){
            try{
                //if(includeCztPcbAbsorption){

                    //effectiveArea = (effarea[pixquad]).get_effective_area((unsigned char) pixx, (unsigned char) pixy,
                    //    expTable.energies);
                //} else {
                    effectiveArea = (effarea[pixquad]).get_effective_area((unsigned char) pixx, (unsigned char) pixy,
                            tempEnergies);
                //}
            } catch(string s){
                LOG(ERROR) << s;
                return EXIT_FAILURE;
            }
            for(ienergy=0; ienergy<expTable.energies.size(); ienergy++){

                //D=sum(wi*ai*fracexp(i))/(sum(ai*fracxp(i)))
                renormOffsets[ienergy] += expTable.weightsArray[i][ienergy]*effectiveArea[ienergy]*exposurefrac_locxy[expTable.locX[i]][expTable.locY[i]];  

                //Take the sum of areas of the pixels which are having non-zero weights
                if(expTable.weightsArray[i][ienergy]!=0)
                    areas[ienergy] += effectiveArea[ienergy]*exposurefrac_locxy[expTable.locX[i]][expTable.locY[i]];

            }

        //fprintf(f1,"%d\t%f\t%f\t%lf\n",i,expTable.weightsArray[i][50],effectiveArea[50],exposurefrac_locxy[expTable.locX[i]][expTable.locY[i]]);

        }
    }
   
     
    for(ienergy=0; ienergy<expTable.energies.size(); ienergy++){
        renormOffsets[ienergy] = renormOffsets[ienergy]/areas[ienergy];
    }
    
    for(i=0; i<expTable.get_nrows(); i++){
        for(ienergy=0; ienergy<expTable.energies.size(); ienergy++){
            if((int)expTable.badpixFlag[i] <= badpixThreshold){
                expTable.weightsArray[i][ienergy] -= renormOffsets[ienergy]; 
            }
        }
    }
    
    expTable.renormOffsetArray.clear();
    expTable.areaArray.clear();
    expTable.renormOffsetArray = renormOffsets;
    expTable.areaArray = areas;
    
    return status;
}
float absco_Al(float ekev) {
    float abscoAl = 0.0;
    float eKev = 0.0;
    double lge, lgabs, x, x2 = 0.0;
    float rhoAl = 2.7;
    
    if(ekev>=1 && ekev<1.56) {
        abscoAl=67.933+1122.14/pow(ekev,3);
        abscoAl=rhoAl*abscoAl;
        return abscoAl;
    }
    
    if(ekev>=1.56 && ekev<2.5){
        x=1.0/ekev;
        x2=x*x;
        abscoAl = -0.9776420 + x * (60.35950 + x2 * (23170.30 - x2 * 21302.5));
        abscoAl = rhoAl*abscoAl;
        return abscoAl;       
    }
    if(ekev>=2.5 && ekev<10) {
        lge = log10 (ekev);
        lgabs = 3.77806+lge*(-0.890608+lge*(-2.29882+lge*0.796674));
        abscoAl = pow(10.0,lgabs);
        abscoAl = rhoAl * abscoAl;
        return abscoAl;       
    }
    if (ekev >= 10 && ekev <= 80) {
        lge = log10 (ekev);
        lgabs = 3.89407 + lge * (-1.29502 + lge * (-1.887 + lge * 0.681882));
        abscoAl = pow(10.0,lgabs);
        abscoAl = rhoAl*abscoAl;
        return abscoAl;
    }
    if (ekev > 80 && ekev <=1000) {
        lge = log10 (ekev);
        lgabs = 6.56821 + lge * (-8.26691 + lge * (3.12323 - lge * 0.410777));
        abscoAl = pow(10.0,lgabs);
        abscoAl = rhoAl*abscoAl;
        return abscoAl;
    }

    return abscoAl;
}

float absco_Ta(float ekev){
    float abscoTa=0.0;
    float eKev=0.0;
    double lge, lgabs=0.0;
    float rhoTa=16.65;

    if (ekev >=1.0 && ekev < 9.88) {
        lge = log10 (ekev);
        lgabs = 4.09765 + lge * (-1.44692 + lge * (-0.700248));
        abscoTa = pow(10.0,lgabs);
        abscoTa = rhoTa * abscoTa;
        return abscoTa;
    }
    if(ekev>=9.88 && ekev<11.14) {
        abscoTa = 754.71 + ekev * (-51.6789);
        abscoTa = rhoTa * abscoTa;        
        return abscoTa;
    }
    if(ekev>=11.14 && ekev<11.68) {
        abscoTa = 801.833 + ekev * (-50.0);
        abscoTa = rhoTa * abscoTa;        
        return abscoTa;
    }
    if(ekev>=11.68 && ekev<67.42) {
        lge = log10 (ekev);
        lgabs = 5.13391 + lge * (-2.5274 + lge * (-0.0274466));
        abscoTa = pow(10.0,lgabs);
        abscoTa = rhoTa*abscoTa;
        return abscoTa;
    }
    if(ekev>=67.42 && ekev<=1000) {
        lge = log10 (ekev);
        lgabs = 2.74145 + lge * (2.52457 + lge * (-2.81744 + lge * 0.513496));
        abscoTa = pow(10.0,lgabs);
        abscoTa = rhoTa * abscoTa;       
        return abscoTa;
    }
    
    return abscoTa;
}

float absco_CZT(float ekev){
    float abscoCZT=0.0;
    float eKev=0.0;
    double lge, lga=0.0;
    float rhoCZT=5.9;   
    
    if(ekev>=5.0 && ekev<26.71) {
        lge = log10 (ekev);
        lga = 4.57449 + lge * (-2.23094 + lge * (-0.208578));
        abscoCZT = pow(10.0,lga);
        abscoCZT = rhoCZT * abscoCZT;
        return abscoCZT;
    }
    if(ekev>=26.71 && ekev<31.81) {
        abscoCZT = 80.7133 - 2.00565 * ekev;
        abscoCZT = rhoCZT * abscoCZT;        
        return abscoCZT;
    }
    if(ekev>=31.81 && ekev<219.65) {
        lge = log10 (ekev);
        lga = -0.332736 + lge * (7.2682 + lge * (-5.55351 + lge * 1.02466));
        abscoCZT = pow(10.0,lga);
        abscoCZT = rhoCZT * abscoCZT;
        return abscoCZT;
    }
    if(ekev>=219.65 && ekev<1000.0) {
        lge = log10 (ekev);
        lga = 26.9552 + lge * (-27.1968 + lge * (8.9179 + lge * (-0.99515)));
        abscoCZT = pow(10.0,lga);
        abscoCZT = rhoCZT * abscoCZT;
        return abscoCZT;
    }
    
    return abscoCZT;
}

float absco_pcb(float ekev){
    float abscoPCB=0.0;
    float eKev=0.0;
    double lge, lga=0.0;
    float rhoPCB=1.412;

    if (ekev >= 0.0 && ekev < 1.739) {
        lga = 3.60207 + lge * (-1.55853 + lge * (-2.03238 + lge * 0.856163));
        abscoPCB = pow(10.0,lga);
        abscoPCB = rhoPCB*abscoPCB;
        return abscoPCB;
    }
    if (ekev >= 1.739 && ekev < 3.5) {
        lga = 0.507307 + lge * (-1.31153 + lge * (0.465186 + lge * (-0.0732322)));
        abscoPCB = pow(10.0,lga);
        abscoPCB = rhoPCB*abscoPCB;
        return abscoPCB;
    }
    if (ekev >= 3.5 && ekev < 5.0) {
        lga = 1.71119 + (-31.6429 / lge)+(72.0848 / (lge * lge));
        abscoPCB = pow(10.0,lga);
        abscoPCB = rhoPCB*abscoPCB;
        return abscoPCB;
    }
    if (ekev >= 5.0) {
        lga = -2.06263 + (6.29368 / lge)+(-23.1632 / (lge * lge));
        abscoPCB = pow(10.0,lga);
        abscoPCB = rhoPCB*abscoPCB;
        return abscoPCB;
    }

    return abscoPCB;
}

float effAbsco(string material, float minEnergy, float maxEnergy, int nIter, int &status){
    int i=0;
    float ekev=minEnergy;
    float effabsco=0.0;
    float delEkev = (maxEnergy-minEnergy)/(nIter-1);
    
    if(maxEnergy>100){
        LOG(ERROR) << "Maximum energy cannot be greater than 100 keV.";
        status = EXIT_FAILURE;
        exit(EXIT_FAILURE);
    }
    if(material=="Al"){
        for(i=0; i<nIter; i++){
            effabsco+=absco_Al(ekev);
            ekev+=delEkev;
        }
        effabsco=effabsco/nIter;
    }
    if(material=="Tan"){
        for(i=0; i<nIter; i++){
            effabsco+=absco_Ta(ekev);
            ekev+=delEkev;
        }
        effabsco=effabsco/nIter;
    }
    if(material=="Czt"){
        for(i=0; i<nIter; i++){
            effabsco+=absco_CZT(ekev);
            ekev+=delEkev;
        }
        effabsco=effabsco/nIter;
    }
    if(material=="Pcb"){
        for(i=0; i<nIter; i++){
            effabsco+=absco_pcb(ekev);
            ekev+=delEkev;
        }
        effabsco=effabsco/nIter;
    }
    
    return effabsco;
}

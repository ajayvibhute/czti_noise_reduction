#include "badpixCALDB.h"


void BadpixTable::pushback_pixFlag(unsigned char detid, unsigned char pixid, 
        unsigned char pixxquad, unsigned char pixyquad, unsigned char locx, unsigned char locy, unsigned char pixflag){
    detID.push_back(detid);
    pixID.push_back(pixid);
    pixXquad.push_back(pixxquad);
    pixYquad.push_back(pixyquad);
    locX.push_back(locx);
    locY.push_back(locy);
    pixFlag.push_back(pixflag);
}

void BadpixTable::reset(){
    detID.clear();
    pixID.clear();
    pixXquad.clear();
    pixYquad.clear();
    locX.clear();
    locY.clear();
    pixFlag.clear();
}

void BadpixTable::resize_table(long size) {
    detID.resize(size, 0);
    pixID.resize(size, 0);
    pixXquad.resize(size, 0);
    pixYquad.resize(size, 0);
    locX.resize(size, 0);
    locY.resize(size, 0);
    pixFlag.resize(size, 0);
}


void BadpixTable::display() {
    map<string,unsigned int>::iterator it = nBadpix.begin();
    int totalPixels=0;
    LOG(INFO) << "Pixel information for quadrant " << quadID;
    DLOG(INFO) << "GOOD:good pixels, SPECBAD: spectroscopically bad pixels";
    DLOG(INFO) << "FLICKERING:flickering pixels, NOISY: noisy pixels";
    DLOG(INFO) << "DEAD: dead/inactive pixels.";
    LOG(INFO) << "-------------------------------------------";
    for(it=nBadpix.begin(); it!=nBadpix.end(); it++){
        LOG(INFO) << it->first << " => " << it->second;
        totalPixels += it->second;
    }
    LOG(INFO) << "TOATAL PIXELS=>" << totalPixels;
    LOG(INFO) << "-------------------------------------------";
    
}

int BadpixTable::generate_badpix_info() {
    int status=0; //status variable
    unsigned int npixels=0;
    map<string, unsigned int>::iterator it=nBadpix.begin();
    
    if(pixFlag.size() < XPIX_QUAD*YPIX_QUAD){
        LOG(WARNING) << "Number of pixels in this badpixel table " << pixFlag.size() <<
                " is less than total number of pixels in a quadrant " << XPIX_QUAD*YPIX_QUAD << ".";
    }
    //Counting good pixel values
    npixels = count(pixFlag.begin(), pixFlag.end(), GOODPIX);
    nBadpix.insert(it, pair<string,unsigned int>("GOOD",npixels));
    //Counting spectroscopically bad pixel values
    npixels = count(pixFlag.begin(), pixFlag.end(), SPECBADPIX);
    nBadpix.insert(it, pair<string,unsigned int>("SPECBAD",npixels));
    //Counting spectroscopically bad pixel values
    npixels = count(pixFlag.begin(), pixFlag.end(), FLICKPIX);
    nBadpix.insert(it, pair<string,unsigned int>("FLICKERING",npixels));
    //Counting spectroscopically bad pixel values
    npixels = count(pixFlag.begin(), pixFlag.end(), NOISYPIX);
    nBadpix.insert(it, pair<string,unsigned int>("NOISY",npixels));
    //Counting spectroscopically bad pixel values
    npixels = count(pixFlag.begin(), pixFlag.end(), DEADPIX);
    nBadpix.insert(it, pair<string,unsigned int>("DEAD",npixels));
    
    return EXIT_SUCCESS;
    
}


int BadpixTable::create_quadrant_badPixMap() {
    int i, j = 0; //counter variables
    int xCoord = 0; // x and y coordinates for badPixMap
    int yCoord = 0;
    vector <unsigned char> vecTemp(NO_PIX_X_PER_DET*NO_DET_X_PER_QUAD, 0);
    vecBadpixMap.clear();
    for (i = 0; i < (NO_PIX_X_PER_DET * NO_DET_X_PER_QUAD); i++) {
        vecBadpixMap.push_back(vecTemp);
    }

    for (i = 0; i < detID.size(); i++) {
        vecBadpixMap[pixYquad[i]][pixXquad[i]] = pixFlag[i];
    }
    return EXIT_SUCCESS;
}

void Badpix::create_badPixMap() {
    int i, j = 0; //counter variables
    int xCoord = 0; // x and y coordinates for badPixMap
    int yCoord = 0;
    int npixels=XPIX_QUAD*YPIX_QUAD;
    string errorMsg="";
    long nrows=0;
    vector <unsigned char> vecTemp(XSIZE, 0);
    badpixMap.clear();
    for (i = 0; i < (YSIZE); i++) {
        badpixMap.push_back(vecTemp);
    }
    //Checking errors in badpix supplied
    if(badpix.size()<4){
        errorMsg = "All 4 quadrants are not present.";
        throw errorMsg;
    }
    for (i = 0; i < badpix.size(); i++) {
        nrows = badpix[i].get_nrows();
        if(nrows!=npixels){
            errorMsg = "Number of rows in quadrant " + itoa(badpix[i].quadID) + " is not equal to 4096.";
            throw errorMsg;
        }
        for (j=0; j< nrows; j++){
            badpixMap[badpix[i].locY[j]][badpix[i].locX[j]] = badpix[i].pixFlag[j];
        }
    }
}

BadpixTable flag_badpix_quadrant(vector<BadpixTable> vecBadpixTable) {
    BadpixTable flaggedTable;
    int nbpixtables=0; //number of badpix tables
    int itable=0;
    int ipixels=0;
    unsigned char pixFlag;
    string errorMsg="";
    nbpixtables = vecBadpixTable.size();
    if(nbpixtables==1){
        return vecBadpixTable[0];
    } else if(nbpixtables < 1){
        errorMsg = "Empty badpix table vector.";
        throw errorMsg;
    }
    
    for(itable=1; itable<nbpixtables; itable++){
        if(vecBadpixTable[0].quadID != vecBadpixTable[itable].quadID){
            errorMsg = "Quadrant IDs are not same for all bad pixel tables.";
            throw errorMsg;
        }
    }
    
    flaggedTable.resize_table(XPIX_QUAD*YPIX_QUAD);
    //Getting pixFlag with highest value.
    for(ipixels=0; ipixels<XPIX_QUAD*YPIX_QUAD; ipixels++){
        pixFlag = vecBadpixTable[0].pixFlag[ipixels];
        for(itable=1; itable<nbpixtables; itable++){
            pixFlag = (pixFlag>vecBadpixTable[itable].pixFlag[ipixels]) ? pixFlag : vecBadpixTable[itable].pixFlag[ipixels];
        }
        flaggedTable.pixFlag[ipixels] = pixFlag;
    }
    flaggedTable.quadID = vecBadpixTable[0].quadID;
    flaggedTable.detID = vecBadpixTable[0].detID;
    flaggedTable.pixID = vecBadpixTable[0].pixID;
    flaggedTable.pixXquad = vecBadpixTable[0].pixXquad;
    flaggedTable.pixYquad = vecBadpixTable[0].pixYquad;
    flaggedTable.locX = vecBadpixTable[0].locX;
    flaggedTable.locY = vecBadpixTable[0].locY;
    
    return flaggedTable;
}

vector <BadpixTable> flag_badpix_files(vector<string> vecBadpixFiles){
    vector <BadpixTable> bptables; //to store badpix table for all 4 quadrants.
    vector <BadpixTable> bptableQ0; //to store badpix table for Quadrant 0.
    vector <BadpixTable> bptableQ1; //to store badpix table for Quadrant 1.
    vector <BadpixTable> bptableQ2; //to store badpix table for Quadrant 2.
    vector <BadpixTable> bptableQ3; //to store badpix table for Quadrant 3.
    Badpix badpix; //to read a badpixel file
    int nfiles=0; //number of badpixfiles.
    int ifile=0; //file index
    int iquad=0; //quadid
    string errorMsg="";
    
    nfiles = vecBadpixFiles.size();
    
    if(nfiles<1){
        errorMsg = "Empty BadpixFiles vector.";
        throw errorMsg;
    }
    else if(nfiles==1){
        if(badpix.read_badpix_file(vecBadpixFiles[0])){
            errorMsg = "Error in reading bad pixel file " + vecBadpixFiles[0];
            throw errorMsg;
        }
        bptables = badpix.get_badpix_tables();
        bptables[0].generate_badpix_info();
        bptables[0].display();
        bptables[1].generate_badpix_info();
        bptables[1].display();
        bptables[2].generate_badpix_info();
        bptables[2].display();
        bptables[3].generate_badpix_info();
        bptables[3].display();
    }
    else{
        for(ifile=0; ifile<nfiles; ifile++) {
            if (badpix.read_badpix_file(vecBadpixFiles[ifile])) {
                errorMsg = "Error in reading bad pixel file " + vecBadpixFiles[ifile];
                throw errorMsg;
            }
            bptables = badpix.get_badpix_tables();
            bptableQ0.push_back(bptables[0]);
            bptableQ1.push_back(bptables[1]);
            bptableQ2.push_back(bptables[2]);
            bptableQ3.push_back(bptables[3]);
        }
        
        bptables.resize(4);
        try {
            bptables[0] = flag_badpix_quadrant(bptableQ0);
            bptables[1] = flag_badpix_quadrant(bptableQ1);
            bptables[2] = flag_badpix_quadrant(bptableQ2);
            bptables[3] = flag_badpix_quadrant(bptableQ3);
        } catch (string s) {
            throw s;
        }
        bptables[0].generate_badpix_info();
        bptables[0].display();
        bptables[1].generate_badpix_info();
        bptables[1].display();
        bptables[2].generate_badpix_info();
        bptables[2].display();
        bptables[3].generate_badpix_info();
        bptables[3].display();
        
    }
    return bptables; 
}
Badpix::Badpix() {

}

int Badpix::read_badpix_file(string badpixelFilename){
    int status=0;
    int i,j=0;
    int iquad=0;
    int npixels=XPIX_QUAD*YPIX_QUAD;
    int colnum=0;
    long nrows=0;
    string errorMsg="";
    char extname[FLEN_VALUE];
    unsigned char quadNo;
    fitsfile *fbadpix; //Pointer to BADPIX CALDB file
    
    
    //resizing badpix vector to store data for all 4 quadrants
    badpix.resize(NUMQUAD);
    
    fits_open_file(&fbadpix, (char*) badpixelFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening badpixel file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    for(iquad=0; iquad<4; iquad++){
       quadToHDU(iquad, extname);
       //reading badpixel file


        fits_movnam_hdu(fbadpix, BINARY_TBL, extname, 0, &status);
        errorMsg = "Error in reading " + (string)extname + " extension in BADPIX file: " + badpixelFilename;
        if (report_error(status, errorMsg)) {
            return EXIT_FAILURE;
        }

        fits_get_num_rows(fbadpix, &nrows, &status);
        errorMsg = "Error in getting number of rows in " + (string)extname + " extension of BADPIX file: " + badpixelFilename;
        if (report_error(status, errorMsg)) {
            return EXIT_FAILURE;
        }  
        
        if(nrows < npixels){
            LOG(ERROR) << "Number of rows are less than " << npixels 
                    << " for extension " << (string)extname;
            return EXIT_FAILURE;
        }
        
        badpix[iquad].reset();
        //DETID column
        if(read_fits_column(fbadpix, "DetID", TBYTE, 1, 1, nrows, badpix[iquad].detID)){
            LOG(ERROR) << "Error in reading DetID column of extension " << extname << 
                    " of badpixel file: " << badpixelFilename;
            return EXIT_FAILURE;
        }
        //PIXID column
        if(read_fits_column(fbadpix, "PixID", TBYTE, 1, 1, nrows, badpix[iquad].pixID)){
            LOG(ERROR) << "Error in reading PixID column of extension " << extname << 
                    " of badpixel file: " << badpixelFilename;
            return EXIT_FAILURE;
        }
        //PIX_FLAG column
        if(read_fits_column(fbadpix, "PIX_FLAG", TBYTE, 1, 1, nrows, badpix[iquad].pixFlag)){
            LOG(ERROR) << "Error in reading PIX_FLAG column of extension " << extname << 
                    " of badpixel file: " << badpixelFilename;
            return EXIT_FAILURE;
        }
        
        badpix[iquad].pixXquad.resize(npixels);
        badpix[iquad].pixYquad.resize(npixels);
        badpix[iquad].locX.resize(npixels);
        badpix[iquad].locY.resize(npixels);
        for(i=0; i<nrows; i++){
            if(generate_locx_locy(badpix[iquad].detID[i], badpix[iquad].pixID[i],
                    (unsigned char) iquad, badpix[iquad].pixXquad[i], badpix[iquad].pixYquad[i],
                    badpix[iquad].locX[i], badpix[iquad].locY[i])){
                LOG(ERROR) << "Error in getting actual bad pixel locations.";
                return EXIT_FAILURE;
            }
        }
        //setting quadrant id for all badpix tables
        badpix[iquad].quadID = iquad;
    }

    badpixMap.clear();
    
    if(read_img(badpixMap, fbadpix, "BADPIX", TBYTE)){
        LOG(WARNING) << "Cannot read BADPIX extension of file: " << badpixelFilename;
    }
    
    
    fits_close_file(fbadpix, &status);
    if (status) {
        LOG(ERROR) << "Error in closing badpixel file: " << badpixelFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}

int Badpix::create_L2badpix_file(string badpixFilename, string badpixTemplate) {
    int status = 0;
    fitsfile *fl2badpix;
    int i = 0, j = 0;
    int iquad = 0;

    if (create_empty_fitsfile(badpixFilename, badpixTemplate)) {
        LOG(ERROR) << "Error in creating empty level-2 event file " << badpixFilename;
        return EXIT_FAILURE;
    }
    
    return (EXIT_SUCCESS);  
}

int Badpix::write_L2badpix_file(string badpixFilename){
    int status=0;
    int iquad=0;
    string errorMsg="";
    fitsfile *fbadpix;
    char extname[FLEN_VALUE];
    
    fits_open_file(&fbadpix, (char*) badpixFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) <<"Error in opening badpixel file " << badpixFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (iquad = 0; iquad < 4; iquad++) {
        quadToHDU(iquad, extname);
        
        fits_movnam_hdu(fbadpix, BINARY_TBL, extname, 0, &status);
        errorMsg = "Error in moving to " + (string) extname + " extension in BADPIX file: " + badpixFilename;
        if (report_error(status, errorMsg)) {
            return EXIT_FAILURE;
        }
        
        if(write_fits_column(fbadpix, "DETID", TBYTE, 1, 1, badpix[iquad].detID)){
            LOG(ERROR) << "Error in writing fits DETID column.";
            return EXIT_FAILURE;
        }
        if(write_fits_column(fbadpix, "PIXID", TBYTE, 1, 1, badpix[iquad].pixID)){
            LOG(ERROR) << "Error in writing fits PIXID column.";
            return EXIT_FAILURE;
        }
        if(write_fits_column(fbadpix, "PIXX", TBYTE, 1, 1, badpix[iquad].pixXquad)){
            LOG(ERROR) << "Error in writing fits PIXX column.";
            return EXIT_FAILURE;
        }
        if(write_fits_column(fbadpix, "PIXY", TBYTE, 1, 1, badpix[iquad].pixYquad)){
            LOG(ERROR) << "Error in writing fits PIXY column.";
            return EXIT_FAILURE;
        }
        if (write_fits_column(fbadpix, "PIX_FLAG", TBYTE, 1, 1, badpix[iquad].pixFlag)) {
            LOG(ERROR) << "Error in writing fits DETID column.";
            return EXIT_FAILURE;
        }
    }

    //creating image for all 4 quadrants and showing error if image cannot be created.
    try{
        create_badPixMap();
    } catch(string s){
        LOG(ERROR) << "Error in creating badpixmap for all quadrants";
        LOG(ERROR) << s;
        return EXIT_FAILURE;
    }
    
    if(write_img(badpixMap, fbadpix, "BADPIX", BYTE_IMG, TBYTE)){
       LOG(ERROR) << "Error in writing image extension.";
    }
    
    
    fits_close_file(fbadpix, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file: " << badpixFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}

int Badpix::find_dead_noisy_pix(string dpiFilename, vector<vector<float> > emap) {
    int status=0;
    
    
    
}

int Badpix::set_badpix_tables(vector<BadpixTable> bptables){
    int ntables=0;
    int itable=0;
    int nQ0=0, nQ1=0, nQ2=0, nQ3=0;
    string errorMsg="";
    vector <int> quadIDs; //quadrant IDs present in bptables
    
    ntables= bptables.size();
    
    if(ntables!=4){
        errorMsg = "Number of badpixel tables expected: 4 found: " + ntables;
        throw errorMsg;
    }
    
    for(itable=0; itable<ntables; itable++){
        quadIDs.push_back(bptables[itable].quadID);
    }

    nQ0 = count(quadIDs.begin(), quadIDs.end(), 0);
    nQ1 = count(quadIDs.begin(), quadIDs.end(), 1);
    nQ2 = count(quadIDs.begin(), quadIDs.end(), 2);
    nQ3 = count(quadIDs.begin(), quadIDs.end(), 3);
    if (nQ0 != 1 || nQ1 != 1 || nQ2 != 1 || nQ3 != 1) {
        errorMsg = "All quadrants are not present in bad pixel tables.";
        errorMsg += "\n                           nQ0=" + itoa(nQ0);
        errorMsg += "\n                           nQ1=" + itoa(nQ1);
        errorMsg += "\n                           nQ2=" + itoa(nQ2);
        errorMsg += "\n                           nQ3=" + itoa(nQ3);
        throw errorMsg;
    }
    
    badpix = bptables;
    
    return EXIT_SUCCESS;
}










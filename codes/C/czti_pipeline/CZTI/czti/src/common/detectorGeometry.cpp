#include "detectorGeometry.h"

CztDetectorGeometry::CztDetectorGeometry() {

}

int CztDetectorGeometry::init_detectors(){
    int status=0;
    int i,j=0;
    
    //initializing vectors xdet0 & ydet0
    xdet0.clear();
    ydet0.clear();
    for(i=0; i<8; i++){
        xdet0.push_back(0.0);
        ydet0.push_back(0.0);
    }
    for(i=0; i<4; i++){
        xdet0[i]=XDMIN+ i*XPITCH;
        ydet0[i]=YDMIN+ i*YPITCH;
        xdet0[7-i]= -xdet0[i] - XLDET;
        ydet0[7-i] = -ydet0[i] - YLDET;
    }
    
    areaGeom = 64*XLDET*YLDET;
    
    return (EXIT_SUCCESS);
}

void CztDetectorGeometry::display_detector_coordinates(){
    int i,j=0; //counter variables
    
    LOG(INFO) << "Displaying detector coordinates [x,y]";
    for(int i=7; i>=0; i--){
        for(int j=0; j<=7; j++){
            cout << "[" << setw(4)<<xdet0[j] << "," << setw(4)<< ydet0[i] << "]";
            cout << "\t";       
        }
        cout << "\n";
    }
    
    LOG(INFO) << "Geometric Area : " << areaGeom;    
}

int CztDetectorGeometry::define_pixel_coordinates(float xdet, float ydet){
    int ipixx, ipixy=0; //counter variables
    int i,j=0; //counter variables
    int status=0; //status variable
    
    //Clearing vectors
    pixx0.clear();
    pixy0.clear();
    pixxw.clear();
    pixyw.clear();
    
    //Initializing vectors pixx0 & pixy0
    for(i=0; i<16; i++){
        pixx0.push_back(0.0);
        pixy0.push_back(0.0);
        pixxw.push_back(0.0);
        pixyw.push_back(0.0);
    }
    
    //Calculating X coordinates for 16 pixel location with corresponding pixel
    //widths in x direction
    for(ipixx=0; ipixx<16; ipixx++){
        pixx0[ipixx] = ipixx*BIG_PIXEL_WIDTH + xdet - 0.015;
        if(ipixx==0){
            pixx0[ipixx] = pixx0[ipixx]+0.015;
        }
        //setting x width of boundary pixels to 0.231cm
        if (ipixx == 0 || ipixx == 15) {
            pixxw[ipixx] = SMALL_PIXEL_WIDTH;
        } else {
            pixxw[ipixx] = BIG_PIXEL_WIDTH;
        }
        
    }

    //Calculating Y coordinates for 16 pixel location with corresponding pixel
    //widths in y direction
    for(ipixy=0; ipixy<16; ipixy++){
        pixy0[ipixy] = ipixy*BIG_PIXEL_WIDTH + ydet - 0.015;
        if(ipixy==0){
            pixy0[ipixy] = pixy0[ipixy]+0.015;
        }
        //setting y width of boundary pixels to 0.231cm
        if (ipixy == 0 || ipixy == 15) {
            pixyw[ipixy] = SMALL_PIXEL_WIDTH;
        } else {
            pixyw[ipixy] = BIG_PIXEL_WIDTH;
        }
    }
    
    return (EXIT_SUCCESS);
}


int CztDetectorGeometry::get_area_map(){
    int status=0;
    int ipixx=0, ipixy=0; //counter variables
    float xdet=0.0;
    float ydet=0.0;
    vector <float> pixWidth;
    vector <float> pixHeight;
    vector <float> tempFullArea;
    // Initializing vectors
    for(ipixx=0; ipixx<XSIZE; ipixx++){
        tempFullArea.push_back(0.0);
        pixWidth.push_back(0.0);
        pixHeight.push_back(0.0);
    }
    areaFullDetector2D.clear();
    for(ipixy=0; ipixy<YSIZE; ipixy++){
        areaFullDetector2D.push_back(tempFullArea);
    }
    //Initialization completed

    for (ipixy = 0; ipixy < YSIZE; ipixy++) {
        for (ipixx = 0; ipixx < XSIZE; ipixx++) {
            switch (ipixy) {
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
                case 127: pixHeight[ipixy] = SMALL_PIXEL_WIDTH;
                    break;
                default: pixHeight[ipixy] = BIG_PIXEL_WIDTH;
            }
            switch (ipixx) {
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
                case 127: pixWidth[ipixx] = SMALL_PIXEL_WIDTH;
                    break;
                default: pixWidth[ipixx] = BIG_PIXEL_WIDTH;
            }
            areaFullDetector2D[ipixy][ipixx] = pixHeight[ipixy]*pixWidth[ipixx];
        }
    }
    return EXIT_SUCCESS;
}
void CztDetectorGeometry::display_pixel_information(string infoType, float xdet, float ydet){
    int ipixx, ipixy=0; //counter variable
    int status=0; //status variable
    float area=0.0; //to store area of each pixel;
    define_pixel_coordinates(xdet, ydet);
    
    if(infoType=="area"){
        LOG(INFO) << "Displaying area of each pixel";
        for(ipixy=15; ipixy>=0; ipixy--){
            for(ipixx=0; ipixx<=15; ipixx++){
                area = pixxw[ipixx]*pixyw[ipixy];
                
                cout<< "[";
                cout<< setprecision(2) << area; 
                cout<< "]\t";
            }
            cout << "\n";
        }
    }
    else if(infoType=="pxloc"){
        LOG(INFO) << "Displaying pixel coordinates";
        for(ipixy=15; ipixy>=0; ipixy--){
            for(ipixx=0; ipixx<=15; ipixx++){
                
                cout<< "[";
                cout<< setw(4) << setprecision(3) << pixx0[ipixx] << "," <<
                        setw(4) << setprecision(3) << pixy0[ipixy]; 
                cout<< "] ";
            }
            cout << "\n";
        }
    }
    
    
    
}



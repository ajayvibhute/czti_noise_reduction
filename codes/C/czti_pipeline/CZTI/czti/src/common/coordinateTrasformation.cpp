#include "coordinateTransformation.h"

int calculate_RA_DEC_Twist(vector <float> vec_nX, vector <float> vec_nY, vector <float> vec_nZ
                            , vector <float> vec_nXt, vector <float> vec_nYt, 
                            vector <float> vec_nZt, float &RA, float &DEC, float &TWIST, int &status){
    float aX, aY, aZ, Tx, Ty, Tz=0.0;
    float tempCT, tempST=0.0;
    
    aX = (accumulate(vec_nX.begin(), vec_nX.end(), 0.0))/ vec_nX.size();
    aY = (accumulate(vec_nY.begin(), vec_nY.end(), 0.0))/vec_nY.size();
    aZ = (accumulate(vec_nZ.begin(), vec_nZ.end(), 0.0))/vec_nZ.size();
    Tx = (accumulate(vec_nXt.begin(), vec_nXt.end(), 0.0))/vec_nXt.size();
    Ty = (accumulate(vec_nYt.begin(), vec_nYt.end(), 0.0))/vec_nYt.size();
    Tz = (accumulate(vec_nZt.begin(), vec_nZt.end(), 0.0))/vec_nZt.size();

    DLOG(INFO) << "AVERAGE OF: 1. nX:" << aX << " 2. nY: " << aY  << " 3. nZ: " << aZ;
    DLOG(INFO) << "AVERAGE OF: 1. nXt:" << Tx << " 2. nYt: " << Ty  << " 3. nZt: " << Tz;
    
    // calculating DEC
    DEC = asin(aZ);
    
    // calculating RA
    if( atan2(aY, aX) < 0) {
        RA = atan2(aY, aX) + 2*M_PI;
    }
    else {
        RA = atan2(aY, aX);
    }
    
    // calculating TWIST angle
    tempCT = Tz*cos(DEC) - Tx*sin(DEC)*cos(RA) - Ty*sin(DEC)*sin(RA);
    tempST = Ty*cos(RA) - Tx*sin(RA);
    
    TWIST = atan2(tempST, tempCT);
    
    DLOG(INFO) << "RA: " << RA << " DEC: " << DEC << " TWIST: " << TWIST;

    return status;
}

int calculate_RA_DEC_Twist(float nX, float nY, float nZ, float nXt, float nYt, float nZt,
                            float &RA, float &DEC, float &TWIST, int &status){
    float tempCT=0.0, tempST=0.0;
    // calculating DEC
    DEC = asin(nZ);

    // calculating RA
    if (atan2(nY, nX) < 0) {
        RA = atan2(nY, nX) + 2 * M_PI;
    } else {
        RA = atan2(nY, nX);
    }

    // calculating TWIST angle
    tempCT = nZt * cos(DEC) - nXt * sin(DEC) * cos(RA) - nYt * sin(DEC) * sin(RA);
    tempST = nYt * cos(RA) - nXt * sin(RA);

    TWIST = atan2(tempST, tempCT);

    DLOG(INFO) << "RA: " << RA << " DEC: " << DEC << " TWIST: " << TWIST;

    return status;
}


int get_inertial_vector(Mkf &mkf, Teldef teldef, double time, float detX, float detY, float detZ, float &nX, float &nY, float &nZ){
    int status=0;
    long startIndex;
    int i,j = 0; //counter variables
    float satX, satY, satZ, icsX, icsY, icsZ=0.0;
    float yawRA, yawDEC, rollRA, rollDEC, pitchRA, pitchDEC =0.0;
    float yawL, yawM, yawN, rollL,rollM, rollN, pitchL, pitchM, pitchN=0.0;
    float alignM11, alignM12, alignM13, alignM21, alignM22, alignM23, alignM31, alignM32, alignM33=0.0;
    float norm=0.0;
    
    startIndex=mkf.get_last_index();
    if (mkf.get_inst_rpy(time, rollRA, rollDEC, pitchRA, pitchDEC, yawRA, yawDEC, startIndex)){
        LOG(ERROR) << "Error in getting roll, pitch & yaw RA & DEC values.";
        return (EXIT_FAILURE);
    }
//    DLOG(INFO) << "rollRA: " << rollRA;
//    DLOG(INFO) << "rollDEC: " << rollDEC;
//    DLOG(INFO) << "pitchRA: " << pitchRA;
//    DLOG(INFO) << "pitchDEC: " << pitchDEC;
//    DLOG(INFO) << "yawRA: " << yawRA;
//    DLOG(INFO) << "yawDEC: " << yawDEC;
    
    // reading alignment values from teldef file
    alignM11 = teldef.get_alignM11();
    alignM12 = teldef.get_alignM12();
    alignM13 = teldef.get_alignM13();
    alignM21 = teldef.get_alignM21();
    alignM22 = teldef.get_alignM22();
    alignM23 = teldef.get_alignM23();
    alignM31 = teldef.get_alignM31();
    alignM32 = teldef.get_alignM32();
    alignM33 = teldef.get_alignM33();
    
    //transforming vectors with components (detX, detY and detZ) in CZTI body frame
    // to satellite reference frame
    satX = alignM11*detX + alignM12*detY + alignM13*detZ;
    satY = alignM21*detX + alignM22*detY + alignM23*detZ;
    satZ = alignM31*detX + alignM32*detY + alignM33*detZ;
    
//    DLOG(INFO) << "satX:" << satX;
//    DLOG(INFO) << "satY:" << satY;
//    DLOG(INFO) << "satZ:" << satZ;
    
    //calculating Direction Cosines
    yawL=cos(yawRA)*cos(yawDEC);
    yawM=sin(yawRA)*cos(yawDEC);
    yawN=sin(yawDEC);
    rollL=cos(rollRA)*cos(rollDEC);
    //rollM=sin(rollRA)*sin(rollDEC);
    rollM=sin(rollRA)*cos(rollDEC);
    rollN=sin(rollDEC);
    pitchL=cos(pitchRA)*cos(pitchDEC);
    pitchM=sin(pitchRA)*cos(pitchDEC);
    pitchN=sin(pitchDEC);
    
//    DLOG(INFO) << "YAWL:" <<yawL;
//    DLOG(INFO) << "YAWM:" <<yawM;
//    DLOG(INFO) << "YAWN:" <<yawN;
//    DLOG(INFO) << "ROLLL:" <<rollL;
//    DLOG(INFO) << "ROLLM:" <<rollM;
//    DLOG(INFO) << "ROLLN:" <<rollN;
//    DLOG(INFO) << "PITCHL:" <<pitchL;
//    DLOG(INFO) << "PITCHM:" <<pitchM;
//    DLOG(INFO) << "PITCHN:" <<pitchN;

    //transforming SATX, SATY and SATZ to inertial coordinates
    icsX=yawL*satX + rollL*satY + pitchL*satZ;
    icsY=yawM*satX + rollM*satY + pitchM*satZ;
    icsZ=yawN*satX + rollN*satY + pitchN*satZ;
    
    norm=sqrt(icsX*icsX + icsY*icsY + icsZ*icsZ);
    
    nX=icsX/norm;
    nY=icsY/norm;
    nZ=icsZ/norm;    
    return status;
} 

int get_body_vector(Mkf mkf, Teldef teldef, float nX, float nY, float nZ, float& detX, 
        float& detY, float& detZ){
    int status=0;
    int i,j = 0; //counter variables
    float satX, satY, satZ =0.0;
    float yawRA, yawDEC, rollRA, rollDEC, pitchRA, pitchDEC =0.0;
    float yawL, yawM, yawN, rollL,rollM, rollN, pitchL, pitchM, pitchN=0.0;
    float alignM11, alignM12, alignM13, alignM21, alignM22, alignM23, alignM31, alignM32, alignM33=0.0;
    float norm=0.0;

    if(mkf.get_avg_rpy(rollRA, rollDEC, pitchRA, pitchDEC, yawRA, yawDEC)){
        LOG(ERROR) << "Error in getting RA & DEC values for roll, pitch and yaw axis";
        status=EXIT_FAILURE;
        return status;
    }

    DLOG(INFO) << "rollRA: " << rollRA;
    DLOG(INFO) << "rollRA: " << rollRA;
    DLOG(INFO) << "pitchRA: " << pitchRA;
    DLOG(INFO) << "pitchDEC: " << pitchDEC;
    DLOG(INFO) << "yawRA: " << yawRA;
    DLOG(INFO) << "yawDEC: " << yawDEC;
    
    //calculating Direction Cosines
    yawL=cos(yawRA)*cos(yawDEC);
    yawM=sin(yawRA)*cos(yawDEC);
    yawN=sin(yawDEC);
    rollL=cos(rollRA)*cos(rollDEC);
    rollM=sin(rollRA)*cos(rollDEC);
    rollN=sin(rollDEC);
    pitchL=cos(pitchRA)*cos(pitchDEC);
    pitchM=sin(pitchRA)*cos(pitchDEC);
    pitchN=sin(pitchDEC);

    // reading alignment values from teldef file
    alignM11 = teldef.get_alignM11();
    alignM12 = teldef.get_alignM12();
    alignM13 = teldef.get_alignM13();
    alignM21 = teldef.get_alignM21();
    alignM22 = teldef.get_alignM22();
    alignM23 = teldef.get_alignM23();
    alignM31 = teldef.get_alignM31();
    alignM32 = teldef.get_alignM32();
    alignM33 = teldef.get_alignM33();

    //transforming from inertial to satellite frame
    satX = yawL*nX + yawM*nY + yawN*nZ;
    satY = rollL*nX + rollM*nY + rollN*nZ;
    satZ = pitchL*nX + pitchM*nY + pitchN*nZ;
    
    //transforming from satellite to CZTI body frame
    detX = alignM11* satX + alignM21*satY + alignM31*satZ;
    detY = alignM12* satX + alignM22*satY + alignM32*satZ;
    detZ = alignM13* satX + alignM23*satY + alignM33*satZ;
    
    norm = sqrt(detX*detX + detY*detY + detZ*detZ);
    
    detX = detX/norm;
    detY = detY/norm;
    detZ = detZ/norm;
    
    return status;
}

int to_nX_nY_nZ(float RA, float DEC, float &nX, float &nY, float &nZ){
    int status=0;
    nX = cos(RA)*cos(DEC);
    nY = sin(RA)*cos(DEC);
    nZ = sin(DEC);
    
    return status;
}

int to_RA_DEC(float nX, float nY, float nZ, float &RA, float &DEC){
    int status=0;
    
    DEC = asin(nZ);
    
    // calculating RA
    if( atan2(nY, nX) < 0) {
        RA = atan2(nY, nX) + 2*M_PI;
    }
    else {
        RA = atan2(nY, nX);
    }
    
    return status;
}

int to_thetaX_thetaY(float detX, float detY, float detZ, float &thetaX, float &thetaY){
    int status=0;
    thetaX = atan(detX/detZ);
    thetaY = atan(detY/detZ);
    
    return status;
}


int to_thetaX_thetaY(string aspectFilename, double RAsource, double DECsource, 
        double &thetaX, double &thetaY){
    AspectFileHandler aspect;
    int status=0; //status variable
    aspect.read_aspect_file(aspectFilename);
    double vX, vY, vZ = 0.0;
    double RA, DEC, TWIST=0.0; // RA, DEC & TWIST of pointing vector to be stored in radians
    RA= aspect.get_RA()*M_PI/180;
    DEC=aspect.get_DEC()*M_PI/180;
    TWIST=aspect.get_TWIST()*M_PI/180;
    //DLOG(INFO) << "Read from aspect file RA:" << RA << "  DEC:" << DEC << " TWIST:" << TWIST;
    vZ = cos(DEC)*cos(DECsource)*cos(RAsource-RA) + sin(DEC)*sin(DECsource);
    vX = -cos(DECsource)*sin(RAsource-RA);
    vY = -( sin(DEC)*cos(DECsource)*cos(RAsource-RA) - sin(DECsource)*cos(DEC));
    
    thetaX = atan((vX*cos(TWIST) + vY*sin(TWIST))/vZ);
    thetaY = atan((vY*cos(TWIST) - vX*sin(TWIST))/vZ);
    return status;
}


int to_detX_detY_detZ(float thetaX, float thetaY, float &detX, float &detY, float &detZ){
    int status=0;
    
    detZ = 1/sqrt(1+ tan(thetaX)*tan(thetaX) + tan(thetaY)*tan(thetaY));
    detX = detZ*tan(thetaX);
    detY = detZ*tan(thetaY);
    
    return status;
}



//DEBUG CODE

//DEBUG CODE END

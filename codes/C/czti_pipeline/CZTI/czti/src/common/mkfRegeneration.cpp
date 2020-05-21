#include "mkfRegeneration.h"
#include "coordinateTransformation.h"
#include "Mvector.h"

using namespace std;


//MKF THRESHOLD CLASS

int MKFthreshold::read_mkf_thresholds(string MKFthrFilename) {
    int status = 0;
    ifstream ifile;
    MKFthrStruct tempMKFthr;
    string line;
    string a, b;

    ifile.open((char*) MKFthrFilename.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        LOG(ERROR) << "Error opening MKF thresholds file " << MKFthrFilename;
        return EXIT_FAILURE;
    }

    LOG(INFO) << "Reading MKF threshold file " << MKFthrFilename;
    while (ifile.good() && (getline(ifile, line))) {
        istringstream iss(line);
        iss >> tempMKFthr.paramaeterName >> tempMKFthr.parFlag >> tempMKFthr.minValue
                >> tempMKFthr.maxValue;
                //LOG(INFO) << " " <<  tempMKFthr.paramaeterName << " " <<  tempMKFthr.parFlag << " " <<  tempMKFthr.minValue
                  //      << " " <<  tempMKFthr.maxValue;
 
 	mkfthr.push_back(tempMKFthr);
    }
    ifile.close(); //closing mkf thresholds file
    LOG(INFO) << "MKF threshold file read and closed. ";

    return status;
}

bool MKFthreshold::verify_mkf_record(MKFrecord mkfr) {
    unsigned int i = 0, j = 0;
    bool flag = false;
    for (i = 0; i < get_npar(); i++) {
        if (mkfthr[i].parFlag == 1) {
            if (mkfthr[i].paramaeterName == "TIME") {
                flag = is_in_range(mkfr.time, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Roll_RA") {
                flag = is_in_range(mkfr.rollRA, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Roll_DEC") {
                flag = is_in_range(mkfr.rollDEC, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "PITCH_RA") {
                flag = is_in_range(mkfr.pitchRA, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "PITCH_DEC") {
                flag = is_in_range(mkfr.pitchDEC, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "YAW_RA") {
                flag = is_in_range(mkfr.yawRA, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "YAW_DEC") {
                flag = is_in_range(mkfr.yawDEC, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "POSX") {
                flag = is_in_range(mkfr.posX, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "POSY") {
                flag = is_in_range(mkfr.posY, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "POSZ") {
                flag = is_in_range(mkfr.posZ, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "SUN_ANGLE") {
                flag = is_in_range(mkfr.sunAngle, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "MOON_ANGLE") {
                flag = is_in_range(mkfr.moonAngle, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "ELV") {
                flag = is_in_range(mkfr.elv, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
                        //Edit by Mithun NPS(10/03/2016) (FOR SUNELV)
           if (mkfthr[i].paramaeterName == "SUNELV") {
                flag = is_in_range(mkfr.sunelev, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }

		   		//Edit by Mithun NPS (05/06/2016) (FOR LAT AND LONG)

           if (mkfthr[i].paramaeterName == "EARTHLAT") {
                flag = is_in_range(mkfr.earthLAT, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }

           if (mkfthr[i].paramaeterName == "EARTHLON") {
                flag = is_in_range(mkfr.earthLON, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
				
            if (mkfthr[i].paramaeterName == "ANG_OFFSET") {
                flag = is_in_range(mkfr.angoffset, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
			
            if (mkfthr[i].paramaeterName == "TIME_SINCE_SAA") {
                flag = is_in_range(mkfr.timeSinceSAA, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
			//edit by ajay. Added CPM_COUNT
            if (mkfthr[i].paramaeterName == "CPM_RATE") {
                flag = is_in_range(mkfr.cpmRate, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }


            //quadrant A
            if (mkfthr[i].paramaeterName == "Q1_MODEID") {
                flag = is_in_range(mkfr.q1modeID, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1_CZTCOUNTER") {
                flag = is_in_range(mkfr.q1CZTcounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1_ALPHACOUNTER") {
                flag = is_in_range(mkfr.q1AlphaCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1_POS_5V_Monitor") {
                flag = is_in_range(mkfr.q1Pos5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1Temperature1") {
                flag = is_in_range(mkfr.q1temp, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1_POS_2DOT5V_Monitor") {
                flag = is_in_range(mkfr.q1Pos2_5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1CZTHV_Monitor") {
                flag = is_in_range(mkfr.q1CZTHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1VetoHV_Monitor") {
                flag = is_in_range(mkfr.q1VetoHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1DVDD") {
                flag = is_in_range(mkfr.q1DVDD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1VetoLLD") {
                flag = is_in_range(mkfr.q1VetoLLD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q1_VetoCounter") {
                flag = is_in_range(mkfr.q1VetoCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            //quadrant B
            if (mkfthr[i].paramaeterName == "Q2_MODEID") {
                flag = is_in_range(mkfr.q2modeID, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2_CZTCOUNTER") {
                flag = is_in_range(mkfr.q2CZTcounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2_ALPHACOUNTER") {
                flag = is_in_range(mkfr.q2AlphaCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2_POS_5V_Monitor") {
                flag = is_in_range(mkfr.q2Pos5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2Temperature1") {
                flag = is_in_range(mkfr.q2temp, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2_POS_2DOT5V_Monitor") {
                flag = is_in_range(mkfr.q2Pos2_5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2CZTHV_Monitor") {
                flag = is_in_range(mkfr.q2CZTHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2VetoHV_Monitor") {
                flag = is_in_range(mkfr.q2VetoHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2DVDD") {
                flag = is_in_range(mkfr.q2DVDD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2VetoLLD") {
                flag = is_in_range(mkfr.q2VetoLLD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q2_VetoCounter") {
                flag = is_in_range(mkfr.q2VetoCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            //quadrant C
            if (mkfthr[i].paramaeterName == "Q3_MODEID") {
                flag = is_in_range(mkfr.q3modeID, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3_CZTCOUNTER") {
                flag = is_in_range(mkfr.q3CZTcounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3_ALPHACOUNTER") {
                flag = is_in_range(mkfr.q3AlphaCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3_POS_5V_Monitor") {
                flag = is_in_range(mkfr.q3Pos5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3Temperature1") {
                flag = is_in_range(mkfr.q3temp, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3_POS_2DOT5V_Monitor") {
                flag = is_in_range(mkfr.q3Pos2_5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3CZTHV_Monitor") {
                flag = is_in_range(mkfr.q3CZTHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3VetoHV_Monitor") {
                flag = is_in_range(mkfr.q3VetoHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3DVDD") {
                flag = is_in_range(mkfr.q3DVDD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3VetoLLD") {
                flag = is_in_range(mkfr.q3VetoLLD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q3_VetoCounter") {
                flag = is_in_range(mkfr.q3VetoCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            //quadrant D
            if (mkfthr[i].paramaeterName == "Q4_MODEID") {
                flag = is_in_range(mkfr.q4modeID, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4_CZTCOUNTER") {
                flag = is_in_range(mkfr.q4CZTcounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4_ALPHACOUNTER") {
                flag = is_in_range(mkfr.q4AlphaCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4_POS_5V_Monitor") {
                flag = is_in_range(mkfr.q4Pos5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4Temperature1") {
                flag = is_in_range(mkfr.q4temp, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4_POS_2DOT5V_Monitor") {
                flag = is_in_range(mkfr.q4Pos2_5vMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4CZTHV_Monitor") {
                flag = is_in_range(mkfr.q4CZTHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4VetoHV_Monitor") {
                flag = is_in_range(mkfr.q4VetoHVMonitor, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4DVDD") {
                flag = is_in_range(mkfr.q4DVDD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4VetoLLD") {
                flag = is_in_range(mkfr.q4VetoLLD, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
            if (mkfthr[i].paramaeterName == "Q4_VetoCounter") {
                flag = is_in_range(mkfr.q4VetoCounter, mkfthr[i].minValue, mkfthr[i].maxValue);
                if (!flag) {
                    break;
                }
            }
        } // end of if (parflag eq 1)
 		//	**********Edit by Mithun NPS (05/06/2016)*************
		else if (mkfthr[i].parFlag >= 2)
		{
			
		   int nconditions=mkfthr[i].parFlag;
		   int ll=0;
		   bool multiflag=false;
		   flag=false;  //	

           if (mkfthr[i].paramaeterName == "EARTHLAT") {
				for(ll=i;ll<i+nconditions;ll++){
                multiflag= is_in_range(mkfr.earthLAT, mkfthr[ll].minValue, mkfthr[ll].maxValue);
				if(multiflag) flag=multiflag;
				}

                if (!flag) {
                    break;
                }
				
            }

           if (mkfthr[i].paramaeterName == "EARTHLON") {
                for(ll=i;ll<i+nconditions;ll++){
                multiflag= is_in_range(mkfr.earthLON, mkfthr[ll].minValue, mkfthr[ll].maxValue);
                if(multiflag) flag=multiflag;
				}
		   
                if (!flag) {
                    break;
                }

			}

			i+=nconditions-1;
		}
    
	}
    return flag;
}

void MKFthreshold::display() {
    int i = 0;
    cout << get_npar();
    for (i = 0; i < get_npar(); i++) {
        LOG(INFO) << mkfthr[i].paramaeterName << ","
                << mkfthr[i].parFlag << ","
                << mkfthr[i].minValue << ","
                << mkfthr[i].maxValue;
    }
}

//MKF THRESHOLD END

//MKF
// Mkf class functions to read, store and get mkf information

Mkf::Mkf() {
    mkfFilename = "";
    optionalParRead = false;
    mkfLastIndex = 0;
}

int Mkf::read_mkf_file(string mkfFileName) {
    int status = 0; // status variable
    int i, j = 0; // counter variables
    string errorMsg="";
    fitsfile *fptr; // Pointer to MKF LEVEL1 FILE.

    // clearing Class vectors
    mkf.clear();
    vec_time.clear();
    vec_Qsat.clear();
    vec_rollRA.clear();
    vec_rollDEC.clear();
    vec_pitchRA.clear();
    vec_pitchDEC.clear();
    vec_yawDEC.clear();
    vec_yawRA.clear();
    vec_posX.clear();
    vec_posY.clear();
    vec_posZ.clear();
    vec_velX.clear();
    vec_velY.clear();
    vec_velZ.clear();
    vec_earthLAT.clear();
    vec_earthLON.clear();
    vec_altitude.clear();
    vec_sunAngle.clear();
    vec_moonAngle.clear();
    vec_elv.clear();
    vec_sunelev.clear();
    vec_angoffset.clear();
    vec_timeSinceSAA.clear();
    vec_cpmRate.clear();
    vec_q1Pos5vMonitor.clear();
    vec_q2Pos5vMonitor.clear();
    vec_q3Pos5vMonitor.clear();
    vec_q4Pos5vMonitor.clear();
    vec_q1CZTcounter.clear();
    vec_q2CZTcounter.clear();
    vec_q3CZTcounter.clear();
    vec_q4CZTcounter.clear();
    vec_q1Pos2_5vMonitor.clear();
    vec_q2Pos2_5vMonitor.clear();
    vec_q3Pos2_5vMonitor.clear();
    vec_q4Pos2_5vMonitor.clear();
    vec_q1CZTHVMonitor.clear();
    vec_q2CZTHVMonitor.clear();
    vec_q3CZTHVMonitor.clear();
    vec_q4CZTHVMonitor.clear();
    vec_q1VetoHVMonitor.clear();
    vec_q2VetoHVMonitor.clear();
    vec_q3VetoHVMonitor.clear();
    vec_q4VetoHVMonitor.clear();
    vec_q1DVDD.clear();
    vec_q2DVDD.clear();
    vec_q3DVDD.clear();
    vec_q4DVDD.clear();
    vec_q1VetoLLD.clear();
    vec_q2VetoLLD.clear();
    vec_q3VetoLLD.clear();
    vec_q4VetoLLD.clear();
    vec_q1VetoCounter.clear();
    vec_q2VetoCounter.clear();
    vec_q3VetoCounter.clear();
    vec_q4VetoCounter.clear();
    vec_q1AlphaCounter.clear();
    vec_q2AlphaCounter.clear();
    vec_q3AlphaCounter.clear();
    vec_q4AlphaCounter.clear();
    vec_q1temp.clear();
    vec_q2temp.clear();
    vec_q3temp.clear();
    vec_q4temp.clear();
    vec_q1modeid.clear();
    vec_q2modeid.clear();
    vec_q3modeid.clear();
    vec_q4modeid.clear();

    //reading mkf file
    fits_open_file(&fptr, mkfFileName.c_str(), READONLY, &status);
    errorMsg = "*** Error in opening MKF file: " + mkfFileName + " ***";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    fits_movnam_hdu(fptr, BINARY_TBL, "MKF", 0, &status);
    errorMsg = "Error in reading MKF extension in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    fits_get_num_rows(fptr, &nrowsMkf, &status);
    LOG(INFO) << "Number of rows in mkf file: " << nrowsMkf;
    errorMsg = "Error in getting number of rows in MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    try {
        read_fits_columnN(fptr, "TIME", TDOUBLE, 1, 1, nrowsMkf, vec_time);
        read_fits_array_columnN(fptr, "Q_SAT", TDOUBLE, 1, 1, nrowsMkf, vec_Qsat);
        read_fits_columnN(fptr, "ROLL_RA", TFLOAT, 1, 1, nrowsMkf, vec_rollRA);
        read_fits_columnN(fptr, "ROLL_DEC", TFLOAT, 1, 1, nrowsMkf, vec_rollDEC);
        read_fits_columnN(fptr, "PITCH_RA", TFLOAT, 1, 1, nrowsMkf, vec_pitchRA);
        read_fits_columnN(fptr, "PITCH_DEC", TFLOAT, 1, 1, nrowsMkf, vec_pitchDEC);
        read_fits_columnN(fptr, "YAW_RA", TFLOAT, 1, 1, nrowsMkf, vec_yawRA);
        read_fits_columnN(fptr, "YAW_DEC", TFLOAT, 1, 1, nrowsMkf, vec_yawDEC);
        read_fits_columnN(fptr, "POSX", TFLOAT, 1, 1, nrowsMkf, vec_posX);
        read_fits_columnN(fptr, "POSY", TFLOAT, 1, 1, nrowsMkf, vec_posY);
        read_fits_columnN(fptr, "POSZ", TFLOAT, 1, 1, nrowsMkf, vec_posZ);
        read_fits_columnN(fptr, "VELX", TFLOAT, 1, 1, nrowsMkf, vec_velX);
        read_fits_columnN(fptr, "VELY", TFLOAT, 1, 1, nrowsMkf, vec_velY);
        read_fits_columnN(fptr, "VELZ", TFLOAT, 1, 1, nrowsMkf, vec_velZ);
        read_fits_columnN(fptr, "EARTHLAT", TFLOAT, 1, 1, nrowsMkf, vec_earthLAT);
        read_fits_columnN(fptr, "EARTHLON", TFLOAT, 1, 1, nrowsMkf, vec_earthLON);
        read_fits_columnN(fptr, "ALTITUDE", TFLOAT, 1, 1, nrowsMkf, vec_altitude);
        read_fits_columnN(fptr, "MOON_ANGLE", TFLOAT, 1, 1, nrowsMkf, vec_moonAngle);
        read_fits_columnN(fptr, "SUN_ANGLE", TFLOAT, 1, 1, nrowsMkf, vec_sunAngle);
        read_fits_columnN(fptr, "ELV", TFLOAT, 1, 1, nrowsMkf, vec_elv);
        read_fits_columnN(fptr, "SUNELV", TFLOAT, 1, 1, nrowsMkf, vec_sunelev);
        read_fits_columnN(fptr, "ANG_OFFSET", TFLOAT, 1, 1, nrowsMkf, vec_angoffset);
        read_fits_columnN(fptr, "TIME_SINCE_SAA", TFLOAT, 1, 1, nrowsMkf, vec_timeSinceSAA);
        read_fits_columnN(fptr, "CPM_RATE", TULONG, 1, 1, nrowsMkf, vec_cpmRate);
        read_fits_columnN(fptr, "Q1_POS_5V_MONITOR", TFLOAT, 1, 1, nrowsMkf, vec_q1Pos5vMonitor);
        read_fits_columnN(fptr, "Q2_POS_5V_MONITOR", TFLOAT, 1, 1, nrowsMkf, vec_q2Pos5vMonitor);
        read_fits_columnN(fptr, "Q3_POS_5V_MONITOR", TFLOAT, 1, 1, nrowsMkf, vec_q3Pos5vMonitor);
        read_fits_columnN(fptr, "Q4_POS_5V_MONITOR", TFLOAT, 1, 1, nrowsMkf, vec_q4Pos5vMonitor);
        read_fits_columnN(fptr, "Q1_POS_2DOT5V_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q1Pos2_5vMonitor);
        read_fits_columnN(fptr, "Q2_POS_2DOT5V_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q2Pos2_5vMonitor);
        read_fits_columnN(fptr, "Q3_POS_2DOT5V_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q3Pos2_5vMonitor);
        read_fits_columnN(fptr, "Q4_POS_2DOT5V_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q4Pos2_5vMonitor);
        read_fits_columnN(fptr, "Q1Temperature1", TFLOAT, 1, 1, nrowsMkf, vec_q1temp);
        read_fits_columnN(fptr, "Q2Temperature1", TFLOAT, 1, 1, nrowsMkf, vec_q2temp);
        read_fits_columnN(fptr, "Q3Temperature1", TFLOAT, 1, 1, nrowsMkf, vec_q3temp);
        read_fits_columnN(fptr, "Q4Temperature1", TFLOAT, 1, 1, nrowsMkf, vec_q4temp);
        read_fits_columnN(fptr, "Q1CZTHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q1CZTHVMonitor);
        read_fits_columnN(fptr, "Q2CZTHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q2CZTHVMonitor);
        read_fits_columnN(fptr, "Q3CZTHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q3CZTHVMonitor);
        read_fits_columnN(fptr, "Q4CZTHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q4CZTHVMonitor);
        read_fits_columnN(fptr, "Q1VetoHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q1VetoHVMonitor);
        read_fits_columnN(fptr, "Q2VetoHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q2VetoHVMonitor);
        read_fits_columnN(fptr, "Q3VetoHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q3VetoHVMonitor);
        read_fits_columnN(fptr, "Q4VetoHV_Monitor", TFLOAT, 1, 1, nrowsMkf, vec_q4VetoHVMonitor);
        read_fits_columnN(fptr, "Q1DVDD", TFLOAT, 1, 1, nrowsMkf, vec_q1DVDD);
        read_fits_columnN(fptr, "Q2DVDD", TFLOAT, 1, 1, nrowsMkf, vec_q2DVDD);
        read_fits_columnN(fptr, "Q3DVDD", TFLOAT, 1, 1, nrowsMkf, vec_q3DVDD);
        read_fits_columnN(fptr, "Q4DVDD", TFLOAT, 1, 1, nrowsMkf, vec_q4DVDD);
        read_fits_columnN(fptr, "Q1VetoLLD", TFLOAT, 1, 1, nrowsMkf, vec_q1VetoLLD);
        read_fits_columnN(fptr, "Q2VetoLLD", TFLOAT, 1, 1, nrowsMkf, vec_q2VetoLLD);
        read_fits_columnN(fptr, "Q3VetoLLD", TFLOAT, 1, 1, nrowsMkf, vec_q3VetoLLD);
        read_fits_columnN(fptr, "Q4VetoLLD", TFLOAT, 1, 1, nrowsMkf, vec_q4VetoLLD);
        read_fits_columnN(fptr, "Q1_VetoCounter", TULONG, 1, 1, nrowsMkf, vec_q1VetoCounter);
        read_fits_columnN(fptr, "Q2_VetoCounter", TULONG, 1, 1, nrowsMkf, vec_q2VetoCounter);
        read_fits_columnN(fptr, "Q3_VetoCounter", TULONG, 1, 1, nrowsMkf, vec_q3VetoCounter);
        read_fits_columnN(fptr, "Q4_VetoCounter", TULONG, 1, 1, nrowsMkf, vec_q4VetoCounter);
        read_fits_columnN(fptr, "Q1_CZT_COUNTER", TULONG, 1, 1, nrowsMkf, vec_q1CZTcounter);
        read_fits_columnN(fptr, "Q2_CZT_COUNTER", TULONG, 1, 1, nrowsMkf, vec_q2CZTcounter);
        read_fits_columnN(fptr, "Q3_CZT_COUNTER", TULONG, 1, 1, nrowsMkf, vec_q3CZTcounter);
        read_fits_columnN(fptr, "Q4_CZT_COUNTER", TULONG, 1, 1, nrowsMkf, vec_q4CZTcounter);
        read_fits_columnN(fptr, "Q1_AlphaCounter", TULONG, 1, 1, nrowsMkf, vec_q1AlphaCounter);
        read_fits_columnN(fptr, "Q2_AlphaCounter", TULONG, 1, 1, nrowsMkf, vec_q2AlphaCounter);
        read_fits_columnN(fptr, "Q3_AlphaCounter", TULONG, 1, 1, nrowsMkf, vec_q3AlphaCounter);
        read_fits_columnN(fptr, "Q4_AlphaCounter", TULONG, 1, 1, nrowsMkf, vec_q4AlphaCounter);
    } catch (ErrorHandler errHandler) {
        logError(errHandler);
        return EXIT_FAILURE;
    }

    try {
        read_fits_column_if_availableN(fptr, "Q1_MODEID", TBYTE, 1, 1, nrowsMkf, vec_q1modeid);
        read_fits_column_if_availableN(fptr, "Q2_MODEID", TBYTE, 1, 1, nrowsMkf, vec_q2modeid);
        read_fits_column_if_availableN(fptr, "Q3_MODEID", TBYTE, 1, 1, nrowsMkf, vec_q3modeid);
        read_fits_column_if_availableN(fptr, "Q4_MODEID", TBYTE, 1, 1, nrowsMkf, vec_q4modeid);
        optionalParRead = true;
    } catch (ErrorHandler errHandler) {
        logError(errHandler);
    }

    fits_close_file(fptr, &status);

    return status;
}

int Mkf::read_mkf_time(fitsfile* fptr) {
    int status = 0;
    int i, j = 0; //counter variables

    int colnum = 0;
    long nrows = 0;
    string errorMsg = "";
    char temp_char[MAX_KEYWORD_SIZE];
    char char_mkfFileName[MAX_KEYWORD_SIZE];
    string mkfFileName;
    double* timeArray; //temporary array to store time column of fits file.

    fits_file_name(fptr, char_mkfFileName, &status);
    errorMsg = "Error in getting name of MKF file.";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    mkfFileName = string(char_mkfFileName);
    fits_movnam_hdu(fptr, BINARY_TBL, "MKF", 0, &status);
    errorMsg = "Error in reading MKF extension in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //reading Time
    fits_get_colnum(fptr, CASEINSEN, "TIME", &colnum, &status);
    errorMsg = "Error in getting column number of TIME column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    timeArray = new double[nrows];
    for (i = 0; i < nrows; i++) {
        timeArray[i] = 0.0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, timeArray, NULL, &status);
    errorMsg = "Error in reading TIME column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    vec_time.clear();

    for (i = 0; i < nrows; i++) {
        vec_time.push_back(timeArray[i]);
    }

    return status;
}

int Mkf::read_mkf_rpy(fitsfile* fptr) {
    int status = 0;
    int i, j = 0; //counter variables
    int colnum = 0;
    long nrows = 0;
    string errorMsg = "";
    char temp_char[MAX_KEYWORD_SIZE];
    char char_mkfFileName[MAX_KEYWORD_SIZE];
    string mkfFileName;

    // temporary array to store information from fits file.
    float* rollRAArray;
    float* rollDECArray;
    float* pitchRAArray;
    float* pitchDECArray;
    float* yawDECArray;
    float* yawRAArray;

    fits_file_name(fptr, char_mkfFileName, &status);
    errorMsg = "Error in getting name of MKF file.";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    mkfFileName = string(char_mkfFileName);
    fits_movnam_hdu(fptr, BINARY_TBL, "MKF", 0, &status);
    errorMsg = "Error in reading MKF extension in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //reading Roll RA
    fits_get_colnum(fptr, CASEINSEN, "ROLL_RA", &colnum, &status);
    errorMsg = "Error in getting column number of ROLL_RA column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    rollRAArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        rollRAArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, rollRAArray, NULL, &status);
    errorMsg = "Error in reading ROLL_RA column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //reading roll DEC
    fits_get_colnum(fptr, CASEINSEN, "ROLL_DEC", &colnum, &status);
    errorMsg = "Error in getting column number of ROLL_DEC column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    rollDECArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        rollDECArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, rollDECArray, NULL, &status);
    errorMsg = "Error in reading ROLL_DEC column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    // reading pitch RA
    fits_get_colnum(fptr, CASEINSEN, "PITCH_RA", &colnum, &status);
    errorMsg = "Error in getting column number of PITCH_RA column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    pitchRAArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        pitchRAArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, pitchRAArray, NULL, &status);
    errorMsg = "Error in reading PITCH_RA column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    // reading pitch DEC
    fits_get_colnum(fptr, CASEINSEN, "PITCH_DEC", &colnum, &status);
    errorMsg = "Error in getting column number of PITCH_DEC column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    pitchDECArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        pitchDECArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, pitchDECArray, NULL, &status);
    errorMsg = "Error in reading PITCH_DEC column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //reading yaw RA
    fits_get_colnum(fptr, CASEINSEN, "YAW_RA", &colnum, &status);
    errorMsg = "Error in getting column number of YAW_RA column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    yawRAArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        yawRAArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, yawRAArray, NULL, &status);
    errorMsg = "Error in reading YAW_RA column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //reading yaw DEC
    fits_get_colnum(fptr, CASEINSEN, "YAW_DEC", &colnum, &status);
    errorMsg = "Error in getting column number of YAW_DEC column in MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    yawDECArray = new float[nrows];
    for (i = 0; i < nrows; i++) {
        yawDECArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, yawDECArray, NULL, &status);
    errorMsg = "Error in reading YAW_DEC column of MKF extension of MKF file: " + mkfFileName;
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    //clearing class vectors
    vec_rollRA.clear();
    vec_rollDEC.clear();
    vec_pitchRA.clear();
    vec_pitchDEC.clear();
    vec_yawDEC.clear();
    vec_yawRA.clear();

    // storing array data in class vectors.
    for (i = 0; i < nrows; i++) {
        vec_rollRA.push_back(rollRAArray[i]);
        vec_rollDEC.push_back(rollDECArray[i]);
        vec_pitchRA.push_back(pitchRAArray[i]);
        vec_pitchDEC.push_back(pitchDECArray[i]);
        vec_yawDEC.push_back(yawDECArray[i]);
        vec_yawRA.push_back(yawRAArray[i]);
    }


    return status;
}

int Mkf::get_inst_rpy(double time, float &rollRA, float &rollDEC, float &pitchRA, float &pitchDEC,
        float &yawRA, float &yawDEC, long startIndex) {
    int status = 0;
    long nrows = vec_time.size();
    double temprollRA = 0.0;
    double temprollDEC = 0.0;
    double temppitchRA = 0.0;
    double temppitchDEC = 0.0;
    double tempyawRA = 0.0;
    double tempyawDEC = 0.0;
    int j = 0; //counter variable
    long i = 0; //counter variable
    long indexMin = nrows + 1;
    long indexMax = nrows + 1;

    for (i = startIndex; i < nrows - 1; i++) {
        if (vec_time[i + 1] > time && vec_time[i] <= time) {
            indexMin = i;
            indexMax = i + 1;
            this->mkfLastIndex = indexMin;
            //DLOG(INFO) << "Event Time " << time << " lies between " << indexMin << " and " << indexMax;
            break;
        }
    }

    if ((indexMin == nrows + 1) || (indexMax == nrows + 1)) {
        LOG(ERROR) << "Time " << setprecision(20) << time << " is outside MKF time range.";
        //return (EXIT_FAILURE);
    }

    //rollRA
    if (vec_rollRA[indexMin] == vec_rollRA[indexMax]) {
        temprollRA = vec_rollRA[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_rollRA[indexMin], (double) vec_rollRA[indexMax], temprollRA);
    }

    //rollDEC
    if (vec_rollDEC[indexMin] == vec_rollDEC[indexMax]) {
        temprollDEC = vec_rollDEC[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_rollDEC[indexMin], (double) vec_rollDEC[indexMax], temprollDEC);
    }

    //pitchRA
    if (vec_pitchRA[indexMin] == vec_pitchRA[indexMax]) {
        temppitchRA = vec_pitchRA[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_pitchRA[indexMin], (double) vec_pitchRA[indexMax], temppitchRA);
    }

    //pitchDEC
    if (vec_pitchDEC[indexMin] == vec_pitchDEC[indexMax]) {
        temppitchDEC = vec_pitchDEC[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_pitchDEC[indexMin], (double) vec_pitchDEC[indexMax], temppitchDEC);
    }

    //yawRA
    if (vec_yawRA[indexMin] == vec_yawRA[indexMax]) {
        tempyawRA = vec_yawRA[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_yawRA[indexMin], (double) vec_yawRA[indexMax], tempyawRA);
    }

    //yawDEC
    if (vec_yawDEC[indexMin] == vec_yawDEC[indexMax]) {
        tempyawDEC = vec_yawDEC[indexMin];
    } else {
        find_interpolated_value(time, vec_time[indexMin], vec_time[indexMax], (double) vec_yawDEC[indexMin], (double) vec_yawDEC[indexMax], tempyawDEC);
    }

    rollRA = temprollRA * M_PI / 180;
    rollDEC = temprollDEC * M_PI / 180;
    pitchRA = temppitchRA * M_PI / 180;
    pitchDEC = temppitchDEC * M_PI / 180;
    yawRA = tempyawRA * M_PI / 180;
    yawDEC = tempyawDEC * M_PI / 180;

    //    DLOG(INFO) << setprecision(12) << "time is: " << time << "   " << vec_time[indexMin] << "   " << vec_time[indexMax];
    //    DLOG(INFO) << setprecision(12) << "ROLL angle is: " << time << "   " << vec_rollRA[indexMin] << "   " << vec_rollRA[indexMax];
    //    DLOG(INFO) << "rollRA: " << rollRA;
    //    DLOG(INFO) << "rollDEC: " << rollDEC;
    //    DLOG(INFO) << "pitchRA: " << pitchRA;
    //    DLOG(INFO) << "pitchDEC: " << pitchDEC;
    //    DLOG(INFO) << "yawRA: " << yawRA;
    //    DLOG(INFO) << "yawDEC: " << yawDEC;

    return status;
}

int Mkf::get_avg_rpy(float &avg_rollRA, float &avg_rollDEC, float &avg_pitchRA, float &avg_pitchDEC, float &avg_yawRA,
        float &avg_yawDEC) {
    int status = 0;
    if (vec_rollRA.size() == 0 || vec_pitchRA.size() == 0 || vec_yawRA.size() == 0 || vec_rollDEC.size() == 0 ||
            vec_pitchDEC.size() == 0 || vec_yawDEC.size() == 0) {
        LOG(ERROR) << "Vectors containing rollRA, rollDEC, pitchRA, pitchDEC, yawRA & yawDEC empty().";
        status = EXIT_FAILURE;
    } else {

        //calculating avg rollRA and rollDEC
        avg_rollRA = (accumulate(vec_rollRA.begin(), vec_rollRA.end(), 0.0)) / vec_rollRA.size();
        avg_rollDEC = (accumulate(vec_rollDEC.begin(), vec_rollDEC.end(), 0.0)) / vec_rollDEC.size();

        //calculating avg pitchRA and pitchDEC
        avg_pitchRA = (accumulate(vec_pitchRA.begin(), vec_pitchRA.end(), 0.0)) / vec_pitchRA.size();
        avg_pitchDEC = (accumulate(vec_pitchDEC.begin(), vec_pitchDEC.end(), 0.0)) / vec_pitchDEC.size();

        //calculating avg yaw RA and yaw DEC
        avg_yawRA = (accumulate(vec_yawRA.begin(), vec_yawRA.end(), 0.0)) / vec_yawRA.size();
        avg_yawDEC = (accumulate(vec_yawDEC.begin(), vec_yawDEC.end(), 0.0)) / vec_yawDEC.size();
    }
    return status;
}

int Mkf::get_filtered_gti(string mkfThrFilename, vector<double>& tstart, vector<double> &tstop) {
    int status = 0;
    long imkf = 0; //mkf record index
    unsigned int ipar = 0; //parameter index
    long nrows = 0;
    unsigned int npar = 0;
    bool switcher = true;
    MKFthreshold mkfthreshold;
    MKFrecord mkfr;
    nrows = get_nrows();
    LOG(INFO) << "nrows: " << nrows;
    //Reading mkf threshold file
    if (mkfthreshold.read_mkf_thresholds(mkfThrFilename)) {
        LOG(ERROR) << "Error in reading mkf threshold file " << mkfThrFilename;
        return EXIT_FAILURE;
    }
    //Validating each record of mkf file
    for (imkf = 0; imkf < nrows; imkf++) {
        if (get_mkf_record(imkf, mkfr)) {
            LOG(ERROR) << "Error in getting MKF record number " << imkf << " from MKF table.";
            return EXIT_FAILURE;
        }
        if (mkfthreshold.verify_mkf_record(mkfr) == switcher) {

            if (switcher == true) {
                tstart.push_back(vec_time[imkf]);
                switcher = false;
                continue;
            } else if (switcher == false) {
                tstop.push_back(vec_time[imkf - 1]);
                switcher = true;
                continue;
            }

        }
    }
    //to get last tstop value
    if (switcher == false) {
        tstop.push_back(vec_time[nrows - 1]);
    }

    return EXIT_SUCCESS;
}

int Mkf::get_mkf_record(long record_no, MKFrecord &mkfr) {
    long nrows = get_nrows();
    if (record_no >= 0 || record_no < nrows) {
        mkfr.time = vec_time[record_no];
        mkfr.rollRA = vec_rollRA[record_no];
        mkfr.rollDEC = vec_rollDEC[record_no];
        mkfr.pitchRA = vec_pitchRA[record_no];
        mkfr.pitchDEC = vec_pitchDEC[record_no];
        mkfr.yawDEC = vec_yawDEC[record_no];
        mkfr.yawRA = vec_yawRA[record_no];
        mkfr.posX = vec_posX[record_no];
        mkfr.posY = vec_posY[record_no];
        mkfr.posZ = vec_posZ[record_no];
        mkfr.velX = vec_velX[record_no];
        mkfr.velY = vec_velY[record_no];
        mkfr.velZ = vec_velZ[record_no];
        mkfr.earthLAT = vec_earthLAT[record_no];
        mkfr.earthLON = vec_earthLON[record_no];
        mkfr.altitude = vec_altitude[record_no];
        mkfr.sunAngle = vec_sunAngle[record_no];
        mkfr.moonAngle = vec_moonAngle[record_no];
        mkfr.elv = vec_elv[record_no];
 		mkfr.sunelev=vec_sunelev[record_no];
		mkfr.angoffset=vec_angoffset[record_no];			
        mkfr.timeSinceSAA = vec_timeSinceSAA[record_no];
        mkfr.cpmRate = vec_cpmRate[record_no];
        mkfr.q1modeID=vec_q1modeid[record_no];
        mkfr.q2modeID=vec_q2modeid[record_no];
        mkfr.q3modeID=vec_q3modeid[record_no];
        mkfr.q4modeID=vec_q4modeid[record_no];
        mkfr.q1Pos5vMonitor=vec_q1Pos5vMonitor[record_no];
        mkfr.q2Pos5vMonitor=vec_q2Pos5vMonitor[record_no];
        mkfr.q3Pos5vMonitor=vec_q3Pos5vMonitor[record_no];
        mkfr.q4Pos5vMonitor=vec_q4Pos5vMonitor[record_no];
        mkfr.q1CZTcounter=vec_q1CZTcounter[record_no];
        mkfr.q2CZTcounter=vec_q2CZTcounter[record_no];
        mkfr.q3CZTcounter=vec_q3CZTcounter[record_no];
        mkfr.q4CZTcounter=vec_q4CZTcounter[record_no];
        mkfr.q1Pos2_5vMonitor=vec_q1Pos2_5vMonitor[record_no];
        mkfr.q2Pos2_5vMonitor=vec_q2Pos2_5vMonitor[record_no];
        mkfr.q3Pos2_5vMonitor=vec_q3Pos2_5vMonitor[record_no];
        mkfr.q4Pos2_5vMonitor=vec_q4Pos2_5vMonitor[record_no];
        mkfr.q1CZTHVMonitor=vec_q1CZTHVMonitor[record_no];
        mkfr.q2CZTHVMonitor=vec_q2Pos2_5vMonitor[record_no];
        mkfr.q3CZTHVMonitor=vec_q3Pos2_5vMonitor[record_no];
        mkfr.q4CZTHVMonitor=vec_q4Pos2_5vMonitor[record_no];
        mkfr.q1VetoHVMonitor=vec_q1VetoHVMonitor[record_no];
        mkfr.q2VetoHVMonitor=vec_q2Pos2_5vMonitor[record_no];
        mkfr.q3VetoHVMonitor=vec_q3Pos2_5vMonitor[record_no];
        mkfr.q4VetoHVMonitor=vec_q4Pos2_5vMonitor[record_no];
        mkfr.q1DVDD=vec_q1DVDD[record_no];
        mkfr.q2DVDD=vec_q2DVDD[record_no];
        mkfr.q3DVDD=vec_q3DVDD[record_no];
        mkfr.q4DVDD=vec_q4DVDD[record_no];
        mkfr.q1VetoLLD=vec_q1VetoLLD[record_no];
        mkfr.q2VetoLLD=vec_q2DVDD[record_no];
        mkfr.q3VetoLLD=vec_q3DVDD[record_no];
        mkfr.q4VetoLLD=vec_q4DVDD[record_no];
        mkfr.q1VetoCounter=vec_q1VetoCounter[record_no];
        mkfr.q2VetoCounter=vec_q2VetoCounter[record_no];
        mkfr.q3VetoCounter=vec_q3VetoCounter[record_no];
        mkfr.q4VetoCounter=vec_q4VetoCounter[record_no];
        mkfr.q1AlphaCounter=vec_q1AlphaCounter[record_no];
        mkfr.q2AlphaCounter=vec_q2AlphaCounter[record_no];
        mkfr.q3AlphaCounter=vec_q3AlphaCounter[record_no];
        mkfr.q4AlphaCounter=vec_q4AlphaCounter[record_no];
        mkfr.q1temp=vec_q1temp[record_no];
        mkfr.q2temp=vec_q2temp[record_no];
        mkfr.q3temp=vec_q3temp[record_no];
        mkfr.q4temp=vec_q4temp[record_no]; 

    } else {
        LOG(ERROR) << "Record number " << record_no << " is greater than number of" <<
                " records in mkf file.";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

//SETTERS MKF

int Mkf::set_last_index(long lastIndex) {
    int status = 0;
    mkfLastIndex = lastIndex;
    return status;
}
// MKF VERIFICATION 

int Mkf::check_time_repetition(string outputFilename) {
    vector<double> uniqueTime;
    int status = 0; //status variable
    long i, nrows, nrowsUnique = 0;
    vector <double>::iterator it;
    int countTime = 0;
    uniqueTime = vec_time;


    //number of rows in vector time
    nrows = vec_time.size();

    //getting unique time values
    it = unique(uniqueTime.begin(), uniqueTime.end());
    uniqueTime.resize(distance(uniqueTime.begin(), it));
    nrowsUnique = uniqueTime.size();

    LOG(INFO) << "Time values: " << nrows;
    LOG(INFO) << "Unique Time values: " << nrowsUnique;

    for (i = 0; i < nrowsUnique; i++) {
        countTime = count(vec_time.begin(), vec_time.end(), uniqueTime[i]);
        if (countTime >= 5) {
            LOG(INFO) << setprecision(25) << uniqueTime[i] << " : " << countTime << endl;
        }
    }
    return status;
}
	
	
int Mkf::regenerate_mkf_values(string attFilename, string lbtFilename, string orbFilename,string hdrFilename){
    int status=0,modeId=0;
    double tmp;
    Attitude attFileHandler;
    LBT lbtFileHandler;
    Orbit orbFileHandler;
    double time=0.0;
    double hkvalues[8]={0};
    long nrows=0;
    long rowNo=0;
    vector <double> vec_startTime; //startTime for alt, orb & lbt
    vector <double> vec_endTime; //endTime for alt, orb & lbt
    double minTime=0.0;
    double maxTime=0.0;
    long iatt=0, ilbt=0, iorb=0, imkf=0;
    double sunAngle=0.0, moonAngle=0.0;
    attStruct attinfo;
    lbtStruct lbtinfo;
    orbStruct orbinfo;
    vector <attStruct> vecAttitude;
    double RAPnt=0.0; //RA pointing (degrees)
    double DECPnt=0.0; //DEC pointing (degrees)
    Q quat; //quaternion to store attitude quaternion; q4 is scalar
    double param;
	long intpart;
 	float fractpart;
    long nrowshdr;	
    start_indexQ0=0;
    start_indexQ1=0;
    start_indexQ2=0;
    start_indexQ3=0;
    fitsfile *fptrhdr;
	double *Time;
	int *DataID=NULL,*HKChannelNo=NULL,*ModeID=NULL,*ErrorCount=NULL,*BootPageNo=NULL;
	long *ADCOutput=NULL,*AlphaCount=NULL,*VetoCount=NULL;
	LOG(INFO)<<"Reading hdr file...."<<hdrFilename.c_str();
	nrowshdr=getNumrows(hdrFilename.c_str(),2);

	Time=(double*)malloc(sizeof(double)*nrowshdr);
	DataID=(int*)malloc(sizeof(int)*nrowshdr);
	HKChannelNo=(int*)malloc(sizeof(int)*nrowshdr);
	ADCOutput=(long*)malloc(sizeof(long)*nrowshdr);
	ModeID=(int*)malloc(sizeof(int)*nrowshdr);
	ErrorCount=(int*)malloc(sizeof(int)*nrowshdr);
	BootPageNo=(int*)malloc(sizeof(int)*nrowshdr);
	AlphaCount=(long*)malloc(sizeof(long)*nrowshdr);
	VetoCount=(long*)malloc(sizeof(long)*nrowshdr);

	fits_open_file(&fptrhdr,hdrFilename.c_str(),READONLY,&status);
	fits_movabs_hdu(fptrhdr,2,NULL, &status);
	if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

	//reading perticular columns
	//cout<<"\n reading hdr file....";
	fits_read_col(fptrhdr,TDOUBLE,1,1,1,nrowshdr,NULL,Time,NULL,&status);
	fits_read_col(fptrhdr,TINT,3,1,1,nrowshdr,NULL,DataID,NULL,&status);
	fits_read_col(fptrhdr,TINT,23,1,1,nrowshdr,NULL,HKChannelNo,NULL,&status);
	fits_read_col(fptrhdr,TLONG,24,1,1,nrowshdr,NULL,ADCOutput,NULL,&status);
	fits_read_col(fptrhdr,TINT,4,1,1,nrowshdr,NULL,ModeID,NULL,&status);
	fits_read_col(fptrhdr,TINT,13,1,1,nrowshdr,NULL,ErrorCount,NULL,&status);
	fits_read_col(fptrhdr,TINT,15,1,1,nrowshdr,NULL,BootPageNo,NULL,&status);
	fits_read_col(fptrhdr,TLONG,26,1,1,nrowshdr,NULL,AlphaCount,NULL,&status);
	fits_read_col(fptrhdr,TLONG,27,1,1,nrowshdr,NULL,VetoCount,NULL,&status);
	if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

    //Reading attitude, lbt and orbit files.
    try{
        attFileHandler.read_attitude_file(attFilename);
        lbtFileHandler.read_lbt_file(lbtFilename);
        orbFileHandler.read_orbit_file(orbFilename);
    } catch(ErrorHandler errHandler){
        throw errHandler;
    }
    vecAttitude=attFileHandler.get_att();
    //Attitude, lbt and orbit files read and data stored.
    
    //Calculating intersection GTIs
    vec_startTime.resize(3, 0.0);
    vec_endTime.resize(3, 0.0);
    vec_startTime[0] = attFileHandler.get_minDataTime();
    vec_startTime[1] = orbFileHandler.get_minDataTime();
    vec_startTime[2] = lbtFileHandler.get_minDataTime();
    vec_endTime[0] = attFileHandler.get_maxDataTime();
    vec_endTime[1] = orbFileHandler.get_maxDataTime();
    vec_endTime[2] = lbtFileHandler.get_maxDataTime();
    minTime = *max_element(vec_startTime.begin(), vec_startTime.end());
    maxTime = *min_element(vec_endTime.begin(), vec_endTime.end()); 
    /*Added by Ajay Vibhute, Dec 29 2015
     * To get the start and stop values
     * */
     tstart=minTime+1;
     tstop=maxTime-2;
    //printf("\n\nmaxtime: %lf",maxTime);
    //printf("\nmintime: %lf",minTime);	
     
     long tstarti=minTime;	
     float tstartf=minTime-tstarti;
     

     long tstopi=maxTime;	
     float tstopf=maxTime-tstopi;

    //printf("\ntstart: %lf",tstart);
    //printf("\n tstop: %lf",tstop);

    //printf("\ntstartf: %f",tstartf);
    //printf("\ntstopf: %f",tstopf);

    //Generating mkf values at an interval of .128*8=1.024 seconds
    //generating mkf values at an interval of 1 sec
    //Time & Q_sat values
    RAPnt = attFileHandler.get_RAPnt();
    DECPnt = attFileHandler.get_DECPnt();
 	intpart=minTime;
    fractpart = minTime-intpart;


	//printf("\nfract:%lf",fractpart);
    //printf("\nintpart%ld",intpart);	
    //Added by Shrikant Chaudhari,March 1 2016
 //  time = intpart+1;
	 //  Changed by Ajay Vibhute, June 6 2016
 	time=intpart+1;//+fractpart;
   	//cout<<"Time:"<<time; 

    //time = minTime;
    
    nrowsMkf=(long)(maxTime-minTime-1);
    mkf.resize(nrowsMkf);
	//cout<<"\n nrowsMkf: "<<nrowsMkf;
    for(imkf=0; imkf<nrowsMkf; imkf++){
	
	//cout<<"\n current imkf: "<<imkf;
	//printf("\n\n time : %lf",time);
	//getting orbit parametrs
        try {
            orbinfo = orbFileHandler.get_interpolated_orbit(time);
            lbtinfo = lbtFileHandler.get_lbt(time);
            attinfo = attFileHandler.get_interpolated_attitude(time);
        } catch (ErrorHandler errHandler) {
            throw errHandler;
        }

        mkf[imkf].time = time;
        mkf[imkf].Qsat = attinfo.qSat;
        quat.q1 = attinfo.qSat[0];
        quat.q2 = attinfo.qSat[1];
        quat.q3 = attinfo.qSat[2];
        quat.q4 = attinfo.qSat[3];
        get_rpy_RA_DEC(quat, &mkf[imkf].rollRA, &mkf[imkf].rollDEC,
                &mkf[imkf].pitchRA, &mkf[imkf].pitchDEC,
                &mkf[imkf].yawRA, &mkf[imkf].yawDEC);
//Added By ajay to add rollRot column	
	mkf[imkf].rollRot=attinfo.rollRot;
        
	mkf[imkf].posX = orbinfo.x;
        mkf[imkf].posY = orbinfo.y;
        mkf[imkf].posZ = orbinfo.z;
        mkf[imkf].velX = orbinfo.vX;
        mkf[imkf].velY = orbinfo.vY;
        mkf[imkf].velZ = orbinfo.vZ;
        mkf[imkf].earthLAT = orbinfo.lat;
        mkf[imkf].earthLON = orbinfo.lon;
        mkf[imkf].altitude = orbinfo.altitude;
        mkf[imkf].elv = calculate_elevation(RAPnt, DECPnt, (double) orbinfo.x, (double) orbinfo.y, (double) orbinfo.z);
	//calculate sun elv
	mkf[imkf].sunelev=calculate_elevation(RAsun,DecSun , (double) orbinfo.x, (double) orbinfo.y, (double) orbinfo.z);

        mkf[imkf].angoffset = calculate_angular_offset(RAPnt, DECPnt, (double) mkf[imkf].rollRA, (double) mkf[imkf].rollDEC);
        try{
            calculate_sun_moon_angle(time, mkf[imkf].rollRA, mkf[imkf].rollDEC, &sunAngle, &moonAngle,&RAsun,&DecSun);
        } catch(ErrorHandler errHandler){
            throw errHandler;
        }
        mkf[imkf].sunAngle = sunAngle; //calculate sunangle later
        mkf[imkf].moonAngle = moonAngle; //calculate moon angle later
        mkf[imkf].timeSinceSAA = 0.0; //calculate time since SAA later
        mkf[imkf].cpmRate = lbtinfo.cpmRate;


//	/*Added by Shrikant to get HK feb 18 2016*/

        //Quadrant A lbt data

 	gethk(Time,DataID,HKChannelNo,ModeID,ErrorCount,BootPageNo,ADCOutput,AlphaCount,VetoCount,0,time,hkparam,&error_count,&boot_page_no,&modeid,&alphacounter,&vetocounter,nrowshdr,&start_indexQ0,&start_indexQ1,&start_indexQ2,&start_indexQ3);        
        // mkf[imkf].q1modeID = (unsigned char) modeId;
  //      mkf[imkf].q1modeID = lbtinfo.q1modeID;
	
	mkf[imkf].q1modeID = modeid;
        mkf[imkf].q1Pos5vMonitor = hkparam[1];
        mkf[imkf].q1temp = hkparam[2];
        mkf[imkf].q1Pos2_5vMonitor = hkparam[3];;
        mkf[imkf].q1VetoHVMonitor = hkparam[5];;
        mkf[imkf].q1CZTHVMonitor =hkparam[4];
        mkf[imkf].q1VetoLLD = hkparam[7];;
        mkf[imkf].q1DVDD = hkparam[6];  
	mkf[imkf].q1VetoCounter = vetocounter;
	mkf[imkf].q1AlphaCounter =alphacounter;

/*        mkf[imkf].q1Pos5vMonitor = lbtinfo.q1Pos5VMonitor;
        mkf[imkf].q1temp = lbtinfo.q1Temperature1;
        mkf[imkf].q1Pos2_5vMonitor = lbtinfo.q1Pos2dot5VMonitor;
        mkf[imkf].q1VetoHVMonitor = lbtinfo.q1VetoHVMonitor;
        mkf[imkf].q1CZTHVMonitor = lbtinfo.q1CZTHVMonitor;
        mkf[imkf].q1VetoCounter = lbtinfo.q1VetoCounter;
        mkf[imkf].q1VetoLLD = lbtinfo.q1vetoLLD;
        mkf[imkf].q1DVDD = lbtinfo.q1DVDD;
      
        mkf[imkf].q1VetoCounter = lbtinfo.q1VetoCounter;
*/      mkf[imkf].q1CZTcounter = lbtinfo.q1CZTCounter;
//	mkf[imkf].q1AlphaCounter = lbtinfo.q1AlphaCounter;


        //Quadrant B lbt data
        
gethk(Time,DataID,HKChannelNo,ModeID,ErrorCount,BootPageNo,ADCOutput,AlphaCount,VetoCount,1,time,hkparam,&error_count,&boot_page_no,&modeid,&alphacounter,&vetocounter,nrowshdr,&start_indexQ0,&start_indexQ1,&start_indexQ2,&start_indexQ3);        
        //mkf[imkf].q2modeID = lbtinfo.q2modeID;
	//mkf[imkf].q2modeID = (unsigned char) modeId;
	mkf[imkf].q2modeID = modeid;
        mkf[imkf].q2Pos5vMonitor = hkparam[1];
        mkf[imkf].q2temp = hkparam[2];;
        mkf[imkf].q2Pos2_5vMonitor = hkparam[3];;
        mkf[imkf].q2VetoHVMonitor = hkparam[5];;
        mkf[imkf].q2CZTHVMonitor =hkparam[4];
        mkf[imkf].q2VetoLLD = hkparam[7];;
        mkf[imkf].q2DVDD = hkparam[6];  
	mkf[imkf].q2VetoCounter = vetocounter;
	mkf[imkf].q2AlphaCounter =alphacounter;

/*	 mkf[imkf].q2Pos5vMonitor = lbtinfo.q2Pos5VMonitor;
        mkf[imkf].q2temp = lbtinfo.q2Temperature1;
        mkf[imkf].q2Pos2_5vMonitor = lbtinfo.q2Pos2dot5VMonitor;
        mkf[imkf].q2VetoHVMonitor = lbtinfo.q2VetoHVMonitor;
        mkf[imkf].q2CZTHVMonitor = lbtinfo.q2CZTHVMonitor;
        mkf[imkf].q2VetoCounter = lbtinfo.q2VetoCounter;
        mkf[imkf].q2VetoLLD = lbtinfo.q2vetoLLD;
        mkf[imkf].q2DVDD = lbtinfo.q2DVDD;

*/        mkf[imkf].q2CZTcounter = lbtinfo.q2CZTCounter;
//        mkf[imkf].q2VetoCounter = lbtinfo.q2VetoCounter;
//        mkf[imkf].q2AlphaCounter = lbtinfo.q2AlphaCounter;

        //Quadrant C lbt data
        //
gethk(Time,DataID,HKChannelNo,ModeID,ErrorCount,BootPageNo,ADCOutput,AlphaCount,VetoCount,2,time,hkparam,&error_count,&boot_page_no,&modeid,&alphacounter,&vetocounter,nrowshdr,&start_indexQ0,&start_indexQ1,&start_indexQ2,&start_indexQ3);
       //mkf[imkf].q3modeID = lbtinfo.q3modeID;
       // mkf[imkf].q3modeID = (unsigned char) modeId;

	mkf[imkf].q3modeID = modeid;
        mkf[imkf].q3Pos5vMonitor = hkparam[1];
        mkf[imkf].q3temp = hkparam[2];;
        mkf[imkf].q3Pos2_5vMonitor = hkparam[3];;
        mkf[imkf].q3VetoHVMonitor = hkparam[5];;
        mkf[imkf].q3CZTHVMonitor =hkparam[4];
        mkf[imkf].q3VetoLLD = hkparam[7];;
        mkf[imkf].q3DVDD = hkparam[6];  
	mkf[imkf].q3VetoCounter = vetocounter;
	mkf[imkf].q3AlphaCounter =alphacounter;


/*	mkf[imkf].q3Pos5vMonitor = lbtinfo.q3Pos5VMonitor;
        mkf[imkf].q3temp = lbtinfo.q3Temperature1;
        mkf[imkf].q3Pos2_5vMonitor = lbtinfo.q3Pos2dot5VMonitor;
        mkf[imkf].q3VetoHVMonitor = lbtinfo.q3VetoHVMonitor;
        mkf[imkf].q3CZTHVMonitor = lbtinfo.q3CZTHVMonitor;
        mkf[imkf].q3VetoCounter = lbtinfo.q3VetoCounter;
        mkf[imkf].q3VetoLLD = lbtinfo.q3vetoLLD;
        mkf[imkf].q2DVDD = lbtinfo.q2DVDD;


        mkf[imkf].q3VetoCounter = lbtinfo.q3VetoCounter;
*/        mkf[imkf].q3CZTcounter = lbtinfo.q3CZTCounter;
//        mkf[imkf].q3AlphaCounter = lbtinfo.q3AlphaCounter;
       //Quadrant D lbt data
gethk(Time,DataID,HKChannelNo,ModeID,ErrorCount,BootPageNo,ADCOutput,AlphaCount,VetoCount,3,time,hkparam,&error_count,&boot_page_no,&modeid,&alphacounter,&vetocounter,nrowshdr,&start_indexQ0,&start_indexQ1,&start_indexQ2,&start_indexQ3);        //mkf[imkf].q3modeID = lbtinfo.q4modeID;
        //mkf[imkf].q4modeID = (unsigned char) modeId;
	
/*	cout<<"\n##### HK parameters ######";
	cout<<"\n Supply Voltage: "<<hkparam[0];
	cout<<"\n Temperature:    "<<hkparam[1];
	cout<<"\n VCCA:		  "<<hkparam[2];
	cout<<"\n CZT HV:	  "<<hkparam[3];
	cout<<"\n VETO HV:	  "<<hkparam[4];
	cout<<"\n DVDD:		  "<<hkparam[5];
	cout<<"\n VETO LLD:	  "<<hkparam[6];
	printf("\n alphacounter:  %ld ",alphacounter);
	printf("\n vetocounter:  %ld ",vetocounter);
*/
	mkf[imkf].q4modeID = modeid;
        mkf[imkf].q4Pos5vMonitor = hkparam[1];
        mkf[imkf].q4temp = hkparam[2];;
        mkf[imkf].q4Pos2_5vMonitor = hkparam[3];;
        mkf[imkf].q4VetoHVMonitor = hkparam[5];;
        mkf[imkf].q4CZTHVMonitor =hkparam[4];
        mkf[imkf].q4VetoLLD = hkparam[7];;
        mkf[imkf].q4DVDD = hkparam[6];  
	mkf[imkf].q4VetoCounter = vetocounter;
	mkf[imkf].q4AlphaCounter =alphacounter;

/*        mkf[imkf].q4modeID = lbtinfo.q4modeID;
        mkf[imkf].q4Pos5vMonitor = lbtinfo.q4Pos5VMonitor;
        mkf[imkf].q4temp = lbtinfo.q4Temperature1;
        mkf[imkf].q4Pos2_5vMonitor = lbtinfo.q4Pos2dot5VMonitor;
        mkf[imkf].q4VetoHVMonitor = lbtinfo.q4VetoHVMonitor;
        mkf[imkf].q4CZTHVMonitor = lbtinfo.q4CZTHVMonitor;
        mkf[imkf].q4VetoLLD = lbtinfo.q4vetoLLD;
        mkf[imkf].q4DVDD = lbtinfo.q4DVDD;

        mkf[imkf].q4VetoCounter = lbtinfo.q4VetoCounter;
*/        mkf[imkf].q4CZTcounter = lbtinfo.q4CZTCounter;
//        mkf[imkf].q4AlphaCounter = lbtinfo.q4AlphaCounter;   
       
        time=time+1.0;
    }
delete []Time; 
delete []DataID;
delete []HKChannelNo;
delete []ADCOutput;
delete []ModeID;
delete []ErrorCount;
delete []BootPageNo;
delete []AlphaCount;
delete []VetoCount;

fits_close_file(fptrhdr,&status);



LOG(INFO)<<"Hk parameter extraction from hdr file is COMPLETED....";
    return EXIT_SUCCESS;
}
int Mkf::write_header(string mkffile,string attfile)
{
	fitsfile *mkf,*att;
	cztHeaderParam param;
	int status=0;	
	param.set_default_values();
	fits_open_file(&mkf, (char*) mkffile.c_str(), READWRITE, &status);
    	if (status) {
    	LOG(ERROR)<<"Error while opening MKF file";
	}

	fits_open_file(&att, (char*) attfile.c_str(), READONLY, &status);
    	if (status) {
    	LOG(ERROR)<<"Error while opening ATT file";
	}
	param.readFromHeader(att);
	param.writeToHeader(mkf);
	param.writeTimekey(tstart,tstop,mkf);
	fits_movnam_hdu(mkf, BINARY_TBL, "MKF", NULL, &status);
    	if (status) {
    		LOG(ERROR)<<"Error while moving to HDU MKF";
	}

	param.readFromHeader(att);
	param.writeToHeader(mkf);
	param.writeTimekey(tstart,tstop,mkf);
 	fits_close_file(att, &status);
	if (status) {
    	LOG(ERROR)<<"Error while closing ATT file";
    	}

	//correcting timedel keyword in primary header
	fits_movabs_hdu(mkf,1,NULL,&status);
	if(status) { printf("\nError(%s:%d):Error in moving to Extension in tct ",__FILE__,__LINE__);fits_report_error(stderr,status); return(EXIT_FAILURE);}
	float timedel=1.0;
	fits_update_key(mkf,TFLOAT,"TIMEDEL",&timedel,"integration time",&status);
	if(status) { printf("\nError(%s:%d):",__FILE__,__LINE__);fits_report_error(stderr,status); return(EXIT_FAILURE);}

 	fits_close_file(mkf, &status);
	if (status) {
    	LOG(ERROR)<<"Error while closing MKF file";
    	}

}
int Mkf::write_regenerated_mkffile(string mkffilename, string rMkfTemplate) {
    int status=0;
    int colnum=0;
    long nrows=0;
    fitsfile *fmkf;
    ErrorHandler errHandler;
    
    if(create_empty_fitsfile(mkffilename, rMkfTemplate)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_CREATING_FITS_FROM_TEMPLATE;
        errHandler.errorMsg = "Error in creating mkf file from corresponding template.";
        throw errHandler;
    }
    
    fits_open_file(&fmkf, (char*) mkffilename.c_str(), READWRITE, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error opening mkf file " + mkffilename;
        throw errHandler;
    }
    
    fits_movnam_hdu(fmkf, BINARY_TBL, "MKF", NULL, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to MKF hdu in file: " + mkffilename;
        throw errHandler;
    }

    fits_get_num_rows(fmkf, &nrows, &status);
    if (status) {
        if (status) {
            fits_read_errmsg(errHandler.fitsErrMsg);
            fits_get_errstatus(status, errHandler.fitsErrTxt);
            errHandler.fitsflag = true;
            errHandler.fitsErrorStatus = status;
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in getting number of rows.";
            throw errHandler;
        }
    }
    //Getting vector of columns to be written in mkf file
    get_vectors();
    
    //TIME
    if(write_fits_column(fmkf, "TIME", TDOUBLE, nrows+1, 1, vec_time)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Time column";
        throw errHandler;
    }
    //QSAT
    if(write_fits_array_column(fmkf, "Q_SAT", TDOUBLE, nrows+1, 1, vec_Qsat)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q_SAT column";
        throw errHandler;
    }
    //ROLL_RA
    if (write_fits_column(fmkf, "ROLL_RA", TFLOAT, nrows + 1, 1, vec_rollRA)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ROLL_RA column";
        throw errHandler;
    }
    //ROLL_DEC
    if (write_fits_column(fmkf, "ROLL_DEC", TFLOAT, nrows + 1, 1, vec_rollDEC)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ROLL_DEC column";
        throw errHandler;
    }
    //PITCH_RA
    if (write_fits_column(fmkf, "PITCH_RA", TFLOAT, nrows + 1, 1, vec_pitchRA)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing PITCH_RA column";
        throw errHandler;
    }
    //PITCH_DEC
    if (write_fits_column(fmkf, "PITCH_DEC", TFLOAT, nrows + 1, 1, vec_pitchDEC)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Time column";
        throw errHandler;
    }
    //YAW_RA
    if (write_fits_column(fmkf, "YAW_RA", TFLOAT, nrows + 1, 1, vec_yawRA)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing YAW_RA column";
        throw errHandler;
    }
    //YAW_DEC
    if (write_fits_column(fmkf, "YAW_DEC", TFLOAT, nrows + 1, 1, vec_yawDEC)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Time column";
        throw errHandler;
    }
    //added by ajay vibhute Dec 17, 2015.
    //To add the rollRot column

    //ROLL_ROT

    if (write_fits_column(fmkf, "ROLL_ROT", TFLOAT, nrows + 1, 1, vec_rollRot)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Time column";
        throw errHandler;
    }

    //POSX
    if (write_fits_column(fmkf, "POSX", TFLOAT, nrows + 1, 1, vec_posX)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing POSX column";
        throw errHandler;
    }
    //POSY
    if (write_fits_column(fmkf, "POSY", TFLOAT, nrows + 1, 1, vec_posY)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing POSY column";
        throw errHandler;
    }
    //POSZ
    if (write_fits_column(fmkf, "POSZ", TFLOAT, nrows + 1, 1, vec_posZ)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing POSZ column";
        throw errHandler;
    }
    //VELX
    if (write_fits_column(fmkf, "VELX", TFLOAT, nrows + 1, 1, vec_velX)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing VELX column";
        throw errHandler;
    }
    //VELY
    if (write_fits_column(fmkf, "VELY", TFLOAT, nrows + 1, 1, vec_velY)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing VELY column";
        throw errHandler;
    }
    //VELZ
    if (write_fits_column(fmkf, "VELZ", TFLOAT, nrows + 1, 1, vec_velZ)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing VELZ column";
        throw errHandler;
    }
    //EARTHLAT
    if (write_fits_column(fmkf, "EARTHLAT", TFLOAT, nrows + 1, 1, vec_earthLAT)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing EARTHLAT column";
        throw errHandler;
    }
    //EARTHLON
    if (write_fits_column(fmkf, "EARTHLON", TFLOAT, nrows + 1, 1, vec_earthLON)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing EARTHLON column";
        throw errHandler;
    }
    //ALTITUDE
    if (write_fits_column(fmkf, "ALTITUDE", TFLOAT, nrows + 1, 1, vec_altitude)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ALTITUDE column";
        throw errHandler;
    }
    //SUN_ANGLE
    if (write_fits_column(fmkf, "SUN_ANGLE", TFLOAT, nrows + 1, 1, vec_sunAngle)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing SUN_ANGLE column";
        throw errHandler;
    }
    //MOON_ANGLE
    if (write_fits_column(fmkf, "MOON_ANGLE", TFLOAT, nrows + 1, 1, vec_moonAngle)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing MOON_ANGLE column";
        throw errHandler;
    }
    //ELV
    if (write_fits_column(fmkf, "ELV", TFLOAT, nrows + 1, 1, vec_elv)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ELV column";
        throw errHandler;
    }
    //write sun elv
    if (write_fits_column(fmkf, "SUNELV", TFLOAT, nrows + 1, 1, vec_sunelev)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ELV column";
        throw errHandler;
    }

    //ANG_OFFSET
    if (write_fits_column(fmkf, "ANG_OFFSET", TFLOAT, nrows + 1, 1, vec_angoffset)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing ANG_OFFSET column";
        throw errHandler;
    }
    //TIME_SINCE_SAA
    if (write_fits_column(fmkf, "TIME_SINCE_SAA", TFLOAT, nrows + 1, 1, vec_timeSinceSAA)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing TIME_SINCE_SAA column";
        throw errHandler;
    }
    //CPM_RATE
    if (write_fits_column(fmkf, "CPM_RATE", TULONG, nrows + 1, 1, vec_cpmRate)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing CPM_RATE column";
        throw errHandler;
    }
    //Q1_MODEID
    if (write_fits_column(fmkf, "Q1_MODEID", TBYTE, nrows + 1, 1, vec_q1modeid)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_MODEID column";
        throw errHandler;
    }
    //Q1_POS_5V_MONITOR
    if (write_fits_column(fmkf, "Q1_POS_5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q1Pos5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_POS_5V_MONITOR column";
        throw errHandler;
    }
    //Q1_CZT_COUNTER
    if (write_fits_column(fmkf, "Q1_CZT_COUNTER", TULONG, nrows + 1, 1, vec_q1CZTcounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_CZT_COUNTER column";
        throw errHandler;
    }
    //Q1TEMPERATURE1
    if (write_fits_column(fmkf, "Q1TEMPERATURE1", TFLOAT, nrows + 1, 1, vec_q1temp)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1TEMPERATURE1 column";
        throw errHandler;
    }
    //Q1_POS_2DOT5V_MONITOR
    if (write_fits_column(fmkf, "Q1_POS_2DOT5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q1Pos2_5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_POS_2DOT5V_MONITOR column";
        throw errHandler;
    }
    //Q1VETOHV_MONITOR
    if (write_fits_column(fmkf, "Q1VETOHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q1VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1VETOHV_MONITOR column";
        throw errHandler;
    }
    //Q1CZTHV_MONITOR
    if (write_fits_column(fmkf, "Q1CZTHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q1CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1CZTHV_MONITOR column";
        throw errHandler;
    }
    //Q1_VETOCOUNTER
    if (write_fits_column(fmkf, "Q1_VETOCOUNTER", TULONG, nrows + 1, 1, vec_q1VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_VETOCOUNTER column";
        throw errHandler;
    }
    //Q1VETOLLD
    if (write_fits_column(fmkf, "Q1VETOLLD", TFLOAT, nrows + 1, 1, vec_q1VetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1VETOLLD column";
        throw errHandler;
    }
    //Q1DVDD
    if (write_fits_column(fmkf, "Q1DVDD", TFLOAT, nrows + 1, 1, vec_q1DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1DEVDD column";
        throw errHandler;
    }
    //Q1_ALPHACOUNTER
    if (write_fits_column(fmkf, "Q1_ALPHACOUNTER", TULONG, nrows + 1, 1, vec_q1AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q1_ALPHACOUNTER column";
        throw errHandler;
    }
    //Q2_MODEID
    if (write_fits_column(fmkf, "Q2_MODEID", TBYTE, nrows + 1, 1, vec_q2modeid)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_MODEID column";
        throw errHandler;
    }
    //Q2_POS_5V_MONITOR
    if (write_fits_column(fmkf, "Q2_POS_5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q2Pos5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_POS_5V_MONITOR column";
        throw errHandler;
    }
    //Q2_CZT_COUNTER
    if (write_fits_column(fmkf, "Q2_CZT_COUNTER", TULONG, nrows + 1, 1, vec_q2CZTcounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_CZT_COUNTER column";
        throw errHandler;
    }
    //Q2TEMPERATURE1
    if (write_fits_column(fmkf, "Q2TEMPERATURE1", TFLOAT, nrows + 1, 1, vec_q2temp)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2TEMPERATURE1 column";
        throw errHandler;
    }
    //Q2_POS_2DOT5V_MONITOR
    if (write_fits_column(fmkf, "Q2_POS_2DOT5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q2Pos2_5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_POS_2DOT5V_MONITOR column";
        throw errHandler;
    }
    //Q2VETOHV_MONITOR
    if (write_fits_column(fmkf, "Q2VETOHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q2VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2VETOHV_MONITOR column";
        throw errHandler;
    }
    //Q2CZTHV_MONITOR
    if (write_fits_column(fmkf, "Q2CZTHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q2CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2CZTHV_MONITOR column";
        throw errHandler;
    }
    //Q2_VETOCOUNTER
    if (write_fits_column(fmkf, "Q2_VETOCOUNTER", TULONG, nrows + 1, 1, vec_q2VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_VETOCOUNTER column";
        throw errHandler;
    }
    //Q2VETOLLD
    if (write_fits_column(fmkf, "Q2VETOLLD", TFLOAT, nrows + 1, 1, vec_q2VetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2VETOLLD column";
        throw errHandler;
    }
    //Q2DVDD
    if (write_fits_column(fmkf, "Q2DVDD", TFLOAT, nrows + 1, 1, vec_q2DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2DEVDD column";
        throw errHandler;
    }
    //Q2_ALPHACOUNTER
    if (write_fits_column(fmkf, "Q2_ALPHACOUNTER", TULONG, nrows + 1, 1, vec_q2AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q2_ALPHACOUNTER column";
        throw errHandler;
    }
    //Q3_MODEID
    if (write_fits_column(fmkf, "Q3_MODEID", TBYTE, nrows + 1, 1, vec_q3modeid)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_MODEID column";
        throw errHandler;
    }
    //Q3_POS_5V_MONITOR
    if (write_fits_column(fmkf, "Q3_POS_5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q3Pos5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_POS_5V_MONITOR column";
        throw errHandler;
    }
    //Q3_CZT_COUNTER
    if (write_fits_column(fmkf, "Q3_CZT_COUNTER", TULONG, nrows + 1, 1, vec_q3CZTcounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_CZT_COUNTER column";
        throw errHandler;
    }
    //Q3TEMPERATURE1
    if (write_fits_column(fmkf, "Q3TEMPERATURE1", TFLOAT, nrows + 1, 1, vec_q3temp)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3TEMPERATURE1 column";
        throw errHandler;
    }
    //Q3_POS_2DOT5V_MONITOR
    if (write_fits_column(fmkf, "Q3_POS_2DOT5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q3Pos2_5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_POS_2DOT5V_MONITOR column";
        throw errHandler;
    }
    //Q3VETOHV_MONITOR
    if (write_fits_column(fmkf, "Q3VETOHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q3VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3VETOHV_MONITOR column";
        throw errHandler;
    }
    //Q3CZTHV_MONITOR
    if (write_fits_column(fmkf, "Q3CZTHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q3CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3CZTHV_MONITOR column";
        throw errHandler;
    }
    //Q3_VETOCOUNTER
    if (write_fits_column(fmkf, "Q3_VETOCOUNTER", TULONG, nrows + 1, 1, vec_q3VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_VETOCOUNTER column";
        throw errHandler;
    }
    //Q3VETOLLD
    if (write_fits_column(fmkf, "Q3VETOLLD", TFLOAT, nrows + 1, 1, vec_q3VetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3VETOLLD column";
        throw errHandler;
    }
    //Q3DVDD
    if (write_fits_column(fmkf, "Q3DVDD", TFLOAT, nrows + 1, 1, vec_q3DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3DEVDD column";
        throw errHandler;
    }
    //Q3_ALPHACOUNTER
    if (write_fits_column(fmkf, "Q3_ALPHACOUNTER", TULONG, nrows + 1, 1, vec_q3AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q3_ALPHACOUNTER column";
        throw errHandler;
    }
    //Q4_MODEID
    if (write_fits_column(fmkf, "Q4_MODEID", TBYTE, nrows + 1, 1, vec_q4modeid)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_MODEID column";
        throw errHandler;
    }
    //Q4_POS_5V_MONITOR
    if (write_fits_column(fmkf, "Q4_POS_5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q4Pos5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_POS_5VMONITOR column";
        throw errHandler;
    }
    //Q4_CZT_COUNTER
    if (write_fits_column(fmkf, "Q4_CZT_COUNTER", TULONG, nrows + 1, 1, vec_q4CZTcounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_CZT_COUNTER column";
        throw errHandler;
    }
    //Q4TEMPERATURE1
    if (write_fits_column(fmkf, "Q4TEMPERATURE1", TFLOAT, nrows + 1, 1, vec_q4temp)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4TEMPERATURE1 column";
        throw errHandler;
    }
    //Q4_POS_2DOT5V_MONITOR
    if (write_fits_column(fmkf, "Q4_POS_2DOT5V_MONITOR", TFLOAT, nrows + 1, 1, vec_q4Pos2_5vMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_POS_2DOT5V_MONITOR column";
        throw errHandler;
    }
    //Q4VETOHV_MONITOR
    if (write_fits_column(fmkf, "Q4VETOHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q4VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4VETOHV_MONITOR column";
        throw errHandler;
    }
    //Q4CZTHV_MONITOR
    if (write_fits_column(fmkf, "Q4CZTHV_MONITOR", TFLOAT, nrows + 1, 1, vec_q4CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4CZTHV_MONITOR column";
        throw errHandler;
    }
    //Q4_VETOCOUNTER
    if (write_fits_column(fmkf, "Q4_VETOCOUNTER", TULONG, nrows + 1, 1, vec_q4VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_VETOCOUNTER column";
        throw errHandler;
    }
    //Q4VETOLLD
    if (write_fits_column(fmkf, "Q4VETOLLD", TFLOAT, nrows + 1, 1, vec_q4VetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4VETOLLD column";
        throw errHandler;
    }
    //Q4DVDD
    if (write_fits_column(fmkf, "Q4DVDD", TFLOAT, nrows + 1, 1, vec_q4DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4DEVDD column";
        throw errHandler;
    }
    //Q4_ALPHACOUNTER
    if (write_fits_column(fmkf, "Q4_ALPHACOUNTER", TULONG, nrows + 1, 1, vec_q4AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_WRITING_FITS_COLUMN;
        errHandler.errorMsg = "Error in writing Q4_ALPHACOUNTER column";
        throw errHandler;
    }
    //All columns written
    
    fits_close_file(fmkf, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error closing regenerated mkf file";
        throw errHandler;
    }
    
    LOG(INFO) << mkffilename << " created successfully";
    
    return EXIT_SUCCESS;
}

int Mkf::get_vectors(){
    int status=0;
    long imkf=0;
    
    //Resizing vectors;
    vec_time.resize(nrowsMkf, 0.0);
    vec_Qsat.resize(nrowsMkf);
    vec_rollDEC.resize(nrowsMkf,0.0);
    vec_rollRA.resize(nrowsMkf,0.0);
    vec_pitchRA.resize(nrowsMkf,0.0);
    vec_pitchDEC.resize(nrowsMkf,0.0);
    vec_yawRA.resize(nrowsMkf,0.0);
    vec_yawDEC.resize(nrowsMkf,0.0);
    vec_rollRot.resize(nrowsMkf,0.0);
   //reset sunelv vector
    vec_sunelev.resize(nrowsMkf,0.0);
    vec_posX.resize(nrowsMkf,0.0);
    vec_posY.resize(nrowsMkf,0.0);
    vec_posZ.resize(nrowsMkf,0.0);
    vec_velX.resize(nrowsMkf,0.0);
    vec_velY.resize(nrowsMkf,0.0);
    vec_velZ.resize(nrowsMkf,0.0);
    vec_earthLAT.resize(nrowsMkf,0.0);
    vec_earthLON.resize(nrowsMkf,0.0);
    vec_altitude.resize(nrowsMkf,0.0);
    vec_sunAngle.resize(nrowsMkf,0.0);
    vec_moonAngle.resize(nrowsMkf,0.0);
    vec_elv.resize(nrowsMkf,0.0);
    vec_angoffset.resize(nrowsMkf,0.0);
    vec_timeSinceSAA.resize(nrowsMkf,0.0);
    vec_cpmRate.resize(nrowsMkf,0.0);
    vec_cpmRate.resize(nrowsMkf,0.0);
    vec_q1modeid.resize(nrowsMkf,0);
    vec_q2modeid.resize(nrowsMkf,0);
    vec_q3modeid.resize(nrowsMkf,0);
    vec_q4modeid.resize(nrowsMkf,0);
    vec_q1Pos5vMonitor.resize(nrowsMkf,0.0);
    vec_q2Pos5vMonitor.resize(nrowsMkf,0.0);
    vec_q3Pos5vMonitor.resize(nrowsMkf,0.0);
    vec_q4Pos5vMonitor.resize(nrowsMkf,0.0);
    vec_q1CZTcounter.resize(nrowsMkf,0);
    vec_q2CZTcounter.resize(nrowsMkf,0);
    vec_q3CZTcounter.resize(nrowsMkf,0);
    vec_q4CZTcounter.resize(nrowsMkf,0);
    vec_q1Pos2_5vMonitor.resize(nrowsMkf,0.0);
    vec_q2Pos2_5vMonitor.resize(nrowsMkf,0.0);
    vec_q3Pos2_5vMonitor.resize(nrowsMkf,0.0);
    vec_q4Pos2_5vMonitor.resize(nrowsMkf,0.0);
    vec_q1CZTHVMonitor.resize(nrowsMkf,0.0);
    vec_q2CZTHVMonitor.resize(nrowsMkf,0.0);
    vec_q3CZTHVMonitor.resize(nrowsMkf,0.0);
    vec_q4CZTHVMonitor.resize(nrowsMkf,0.0);
    vec_q1VetoHVMonitor.resize(nrowsMkf,0.0);
    vec_q2VetoHVMonitor.resize(nrowsMkf,0.0);
    vec_q3VetoHVMonitor.resize(nrowsMkf,0.0);
    vec_q4VetoHVMonitor.resize(nrowsMkf,0.0);
    vec_q1DVDD.resize(nrowsMkf,0.0);
    vec_q2DVDD.resize(nrowsMkf,0.0);
    vec_q3DVDD.resize(nrowsMkf,0.0);
    vec_q4DVDD.resize(nrowsMkf,0.0);
    vec_q1VetoLLD.resize(nrowsMkf,0.0);
    vec_q2VetoLLD.resize(nrowsMkf,0.0);
    vec_q3VetoLLD.resize(nrowsMkf,0.0);
    vec_q4VetoLLD.resize(nrowsMkf,0.0);
    vec_q1VetoCounter.resize(nrowsMkf,0);
    vec_q2VetoCounter.resize(nrowsMkf,0);
    vec_q3VetoCounter.resize(nrowsMkf,0);
    vec_q4VetoCounter.resize(nrowsMkf,0);
    vec_q1AlphaCounter.resize(nrowsMkf,0);
    vec_q2AlphaCounter.resize(nrowsMkf,0);
    vec_q3AlphaCounter.resize(nrowsMkf,0);
    vec_q4AlphaCounter.resize(nrowsMkf,0);
    vec_q1temp.resize(nrowsMkf,0.0);
    vec_q2temp.resize(nrowsMkf,0.0);
    vec_q3temp.resize(nrowsMkf,0.0);
    vec_q4temp.resize(nrowsMkf,0.0);
    
    //assigning mkf values
    for(imkf=0; imkf<nrowsMkf; imkf++) {
        vec_time[imkf] = mkf[imkf].time;
        vec_Qsat[imkf] = mkf[imkf].Qsat;
        vec_rollDEC[imkf] = mkf[imkf].rollDEC;
        vec_rollRA[imkf] = mkf[imkf].rollRA;
        vec_pitchRA[imkf] = mkf[imkf].pitchRA;
        vec_pitchDEC[imkf] = mkf[imkf].pitchDEC;
        vec_yawRA[imkf] = mkf[imkf].yawRA;
        vec_yawDEC[imkf] = mkf[imkf].yawDEC;
	vec_rollRot[imkf]=mkf[imkf].rollRot;
   	vec_sunelev[imkf]=mkf[imkf].sunelev;
	//assign the value
        vec_posX[imkf] = mkf[imkf].posX;
        vec_posY[imkf] = mkf[imkf].posY;
        vec_posZ[imkf] = mkf[imkf].posZ;
        vec_velX[imkf] = mkf[imkf].velX;
        vec_velY[imkf] = mkf[imkf].velY;
        vec_velZ[imkf] = mkf[imkf].velZ ;
        vec_earthLAT[imkf] = mkf[imkf].earthLAT ;
        vec_earthLON[imkf] = mkf[imkf].earthLON ;
        vec_altitude[imkf] = mkf[imkf].altitude ;
        vec_sunAngle[imkf] = mkf[imkf].sunAngle ;
        vec_moonAngle[imkf] = mkf[imkf].moonAngle ;
        vec_elv[imkf] = mkf[imkf].elv ;
        vec_angoffset[imkf] = mkf[imkf].angoffset ;
        vec_timeSinceSAA[imkf] = mkf[imkf].timeSinceSAA ;
        vec_cpmRate[imkf] = mkf[imkf].cpmRate ;
        vec_q1modeid[imkf] = mkf[imkf].q1modeID;
        vec_q2modeid[imkf] = mkf[imkf].q2modeID;
        vec_q3modeid[imkf] = mkf[imkf].q3modeID;
        vec_q4modeid[imkf] = mkf[imkf].q4modeID;
        vec_q1Pos5vMonitor[imkf] = mkf[imkf].q1Pos5vMonitor ;
        vec_q2Pos5vMonitor[imkf] = mkf[imkf].q1Pos5vMonitor ;
        vec_q3Pos5vMonitor[imkf] = mkf[imkf].q1Pos5vMonitor ;
        vec_q4Pos5vMonitor[imkf] = mkf[imkf].q1Pos5vMonitor ;
        vec_q1CZTcounter[imkf] = mkf[imkf].q1CZTcounter ;
        vec_q2CZTcounter[imkf] = mkf[imkf].q2CZTcounter ;
        vec_q3CZTcounter[imkf] = mkf[imkf].q3CZTcounter ;
        vec_q4CZTcounter[imkf] = mkf[imkf].q4CZTcounter ;
        vec_q1Pos2_5vMonitor[imkf] = mkf[imkf].q1Pos2_5vMonitor ;
        vec_q2Pos2_5vMonitor[imkf] = mkf[imkf].q2Pos2_5vMonitor ;
        vec_q3Pos2_5vMonitor[imkf] = mkf[imkf].q3Pos2_5vMonitor ;
        vec_q4Pos2_5vMonitor[imkf] = mkf[imkf].q4Pos2_5vMonitor ;
        vec_q1CZTHVMonitor[imkf] = mkf[imkf].q1CZTHVMonitor ;
        vec_q2CZTHVMonitor[imkf] = mkf[imkf].q2CZTHVMonitor ;
        vec_q3CZTHVMonitor[imkf] = mkf[imkf].q3CZTHVMonitor ;
        vec_q4CZTHVMonitor[imkf] = mkf[imkf].q4CZTHVMonitor ;
        vec_q1VetoHVMonitor[imkf] = mkf[imkf].q1VetoHVMonitor ;
        vec_q2VetoHVMonitor[imkf] = mkf[imkf].q2VetoHVMonitor ;
        vec_q3VetoHVMonitor[imkf] = mkf[imkf].q3VetoHVMonitor ;
        vec_q4VetoHVMonitor[imkf] = mkf[imkf].q4VetoHVMonitor ;
        vec_q1DVDD[imkf] = mkf[imkf].q1DVDD ;
        vec_q2DVDD[imkf] = mkf[imkf].q2DVDD ;
        vec_q3DVDD[imkf] = mkf[imkf].q3DVDD ;
        vec_q4DVDD[imkf] = mkf[imkf].q4DVDD ;
        vec_q1VetoLLD[imkf] = mkf[imkf].q1VetoLLD ;
        vec_q2VetoLLD[imkf] = mkf[imkf].q2VetoLLD ;
        vec_q3VetoLLD[imkf] = mkf[imkf].q3VetoLLD ;
        vec_q4VetoLLD[imkf] = mkf[imkf].q4VetoLLD ;
        vec_q1VetoCounter[imkf] = mkf[imkf].q1VetoCounter ;
        vec_q2VetoCounter[imkf] = mkf[imkf].q2VetoCounter ;
        vec_q3VetoCounter[imkf] = mkf[imkf].q3VetoCounter ;
        vec_q4VetoCounter[imkf] = mkf[imkf].q4VetoCounter ;
        vec_q1AlphaCounter[imkf] = mkf[imkf].q1AlphaCounter ;
        vec_q2AlphaCounter[imkf] = mkf[imkf].q2AlphaCounter ;
        vec_q3AlphaCounter[imkf] = mkf[imkf].q3AlphaCounter ;
        vec_q4AlphaCounter[imkf] = mkf[imkf].q4AlphaCounter ;
        vec_q1temp[imkf] = mkf[imkf].q1temp ;
        vec_q2temp[imkf] = mkf[imkf].q2temp ;
        vec_q3temp[imkf] = mkf[imkf].q3temp ;
        vec_q4temp[imkf] = mkf[imkf].q4temp ;
    }
    return EXIT_SUCCESS;
}

double calculate_elevation(double RApnt, double DECpnt, double x, double y, double z) {
    MVector pvec(RApnt, DECpnt); //Pointing vector
    MVector radvec(-x,-y,-z); //vector pointing towards centre of earth
    double thetaS=0.0, thetaE=0.0, elvDeg=0.0; //all angles are in degrees
    thetaS = get_angle_bw_vectors(pvec, radvec);
    thetaE = asin(RADEARTH/radvec.get_magnitude()) * TODEG;
    elvDeg = thetaS - thetaE;
    return elvDeg;  
}

double calculate_angular_offset(double RApnt, double DECpnt, double rollRA, double rollDEC){
    MVector pvec(RApnt, DECpnt); //Pointing vector
    MVector avec(rollRA, rollDEC); //Roll direction vector
    double angularOffsetDeg=0.0; //Angle between pointing and roll direction.
    
    //calculating angular offset (in degrees)
    angularOffsetDeg = get_angle_bw_vectors(pvec, avec);
    return angularOffsetDeg;
}

int sunmoon(double utTime, double *RAsun, double* DecSun, 
        double* RAmoon, double *DecMoon){
    double dayref=(utTime/86400.0)+3654.0;
    double d0=dayref;
    double ws=0.0, es=0.0, Ms=0.0, ob=0.0;
    double Ls=0.0, E0s=0.0;
    double x=0.0, y=0.0, r=0.0, v=0.0;
    double ye=0.0, ze=0.0;
    double lon=0.0, lat=0.0;
    
    double N=0.0, incl=0.0, wm=0.0, am=0.0, em=0.0, Mm=0.0;
    double E0m=0.0, Ediff=0.0;
    double E1m=0.0;
    double Na=0.0, vwa=0.0, ia=0.0, xec=0.0, yec=0.0, zec=0.0; 
    double dr=TORAD;
    double Lm=0.0, D=0.0, F=0.0;
    double dlon=0.0, dlat=0.0, delr=0.0;
    double yeq=0.0, zeq=0.0;
    double fRAdeg=0.0, fDecdeg=0.0;
    //Sun
    ws=282.944+(4.70935e-5)*d0;
    es=0.01609-1.151e-9*d0;
    Ms=356.0470+0.9856002585*d0;
    ob=23.4393-3.63e-7*d0;
    //
    Ls=ws+Ms;
    E0s=Ms+es*sin(dr*Ms)*(1.0+es*cos(dr*Ms))/dr;
    //
    x=cos(dr*E0s)-es;
    y=sin(dr*E0s)*sqrt(1.0-es*es);
    r=sqrt(x*x+y*y);
    if (atan2(y, x) < 0) {
        v = atan2(y, x)/dr + 360;
    } else {
        v = atan2(y,x)/dr;
    }
    //
    lon=v+ws;
    x=r*cos(dr*lon);
    y=r*sin(dr*lon);
    ye=y*cos(dr*ob);
    ze=y*sin(dr*ob);
    r=sqrt(x*x+ye*ye+ze*ze);
    //
    if (atan2(y, x) < 0) {
        *RAsun = atan2(ye, x) / dr + 360;
    } else {
        *RAsun = atan2(ye, x) / dr;
    }
    *DecSun = asin(ze/r)/dr; //in degrees
    
    //Moon
    N=125.1228-0.0529538083*d0;
    incl=5.1454;
    wm=318.0634+0.1643573223*d0;
    am=60.2666;
    em=0.054900;
    Mm=115.3654+13.0649929509*d0;
    //
    E0m=Mm+em*sin(Mm*dr)*(1.0+em*cos(Mm*dr))/dr;
    Ediff=1.0;
    do {
        E1m = E0m - (E0m - em * sin(E0m * dr) / dr - Mm) / (1.0 - em * cos(E0m * dr));
        Ediff = abs(E0m - E1m);
        E0m = E1m;
    } while (Ediff > 0.005);

    x=am*(cos(E0m*dr)-em);
    y=sin(E0m*dr)*am*sqrt(1.0-em*em);
    r=sqrt(x*x+y*y);
    if (atan2(y, x) < 0) {
        v = atan2(y, x) / dr + 360;
    } else {
        v = atan2(y, x) / dr;
    }
    //
    Na=N*dr;
    vwa=(v+wm)*dr;
    ia=incl*dr;
    xec=r*(cos(Na)*cos(vwa)-sin(Na)*sin(vwa)*cos(ia));
    yec=r*(sin(Na)*cos(vwa)+cos(Na)*sin(vwa)*cos(ia));
    zec=r*sin(vwa)*sin(ia);
    if (atan2(yec, xec) < 0) {
        lon = atan2(yec, xec) / dr + 360;
    } else {
        lon = atan2(yec, xec) / dr;
    }
    lat=asin(zec/r)/dr;
    //Moon perturbations
    Lm=N+wm+Mm;
    D=Lm-Ls;
    F=Lm-N;
    //
    dlon=-1.274*sin(dr*(Mm-2.*D))+0.658*sin(2.*D*dr)-0.186*sin(Ms*dr);
    dlon=dlon-0.059*sin(2.*dr*(Mm-D))-0.057*sin(dr*(Mm-2*D+Ms));
    dlon=dlon+0.053*sin(dr*(Mm+2*D))+0.046*sin(dr*(2*D-Ms));
    dlon=dlon+0.041*sin(dr*(Mm-Ms))-0.035*sin(dr*D)-0.031*sin(dr*(Mm+Ms));
    dlon=dlon-0.015*sin(2*dr*(F-D))+0.011*sin(dr*(Mm-4.*D));
    //
    dlat=-0.173*sin(dr*(F-2*D))-0.055*sin(dr*(Mm-F-2*D));
    dlat=dlat-0.046*sin(dr*(Mm+F-2*D))+0.033*sin(dr*(F+2*D));
    dlat = dlat + 0.017 * sin(dr * (2 * Mm + F));
    //
    delr = -0.58 * cos(dr * (Mm - 2 * D)) - 0.46 * cos(2 * dr * D);
    //
    lon = lon + dlon;
    lat = lat + dlat;
    r = r + delr;
    //    
    zec = r * sin(dr * lat);
    yec = r * cos(dr * lat) * sin(dr * lon);
    xec = r * cos(dr * lat) * cos(dr * lon);
    //
    yeq = yec * cos(ob * dr) - zec * sin(ob * dr);
    zeq = yec * sin(ob * dr) + zec * cos(ob * dr);
    if (atan2(yeq, xec) < 0) {
        *RAmoon = atan2(yeq, xec) / dr + 360;
    } else {
        *RAmoon = atan2(yeq, xec) / dr;
    }
    *DecMoon = asin(zeq / r) / dr;

    return EXIT_SUCCESS;
}
int calculate_sun_moon_angle(double utTime, double rollRA, double rollDEC, double* sunAngle, double* moonAngle,double *OutRAsun,double *OutDecSun)
{
    int status=0;
    double dayref = (utTime / 86400.0) + 3654.0;
    double RAmoon=0.0, DecMoon=0.0;
    double fRAdeg=0.0, fDecdeg=0.0;
    double dr=TORAD;
    double RAsun=0.0,DecSun=0.0;
    sunmoon(utTime, &RAsun, &DecSun, &RAmoon, &DecMoon);
    *OutRAsun=RAsun;
    *OutDecSun=DecSun;
    MVector svec(RAsun, DecSun); //Pointing vector to sun
    MVector mvec(RAmoon, DecMoon); //Pointing vector to moon
    
    //PRECESSING rollRA & rollDec from J200 to date
    fRAdeg = 1.14077e-5 * (3.075 + 1.336 * sin(rollRA * dr) * tan(rollDEC * dr));
    fDecdeg = 1.5207e-5 * cos(rollRA * dr);
    rollRA = rollRA + fRAdeg*(dayref-1);
    rollDEC = rollDEC + fDecdeg*(dayref-1);
    
    MVector avec(rollRA, rollDEC); //Roll direction vector

    //calculating sun angle and moon angle (in degrees)
    *sunAngle = get_angle_bw_vectors(svec, avec);
    *moonAngle = get_angle_bw_vectors(mvec, avec);
//    LOG(INFO) << "--------------------------------------";
//    LOG(INFO) << "Sun";
//    svec.display();
//    LOG(INFO) << "Moon";
//    mvec.display();
//    LOG(INFO) << "Precessed Pointing direction (roll)";
//    avec.display();
//    LOG(INFO) << "UTTIME: " << setprecision(20) << utTime;
//    LOG(INFO) << RAsun << " " << DecSun << " " << *sunAngle;
//    LOG(INFO) << RAmoon << " " << DecMoon << " " << *moonAngle;
//    LOG(INFO) << "--------------------------------------";
    
    
    return EXIT_SUCCESS;
}
/*Code extract hk parameters from hdr file
*/
//added by shrikant
int gethk(double Time[],int DataID[],int HKChannelNo[],int ModeID[],int ErrorCount[],int BootPageNo[],long ADCOutput[],long AlphaCount[],long VetoCount[], int quadrant,double search_time,float hkparam[],int *error_count,int *boot_page_no,int *modeid,long *alphacounter,long *vetocounter,long nrowshdr,long *start_indexQ0,long *start_indexQ1,long *start_indexQ2,long *start_indexQ3)
{	
	int j;
	long found_time_quadrant;
	long i,found_time;
	int hkval,flag=1;
	int count=0;
	long adc[8];
	long Qlast=0,Qfirst=0;

adc[0]=-1;
adc[1]=-1;
adc[2]=-1;
adc[3]=-1;
adc[4]=-1;
adc[5]=-1;
adc[6]=-1;
adc[7]=-1;

if(quadrant==0){
i=*start_indexQ0;
}

if(quadrant==1){
i=*start_indexQ1;
}

if(quadrant==2){
i=*start_indexQ2;
}

if(quadrant==3){
i=*start_indexQ3;
}


//getting end time quadrantwise
Qlast=nrowshdr-1;
while(DataID[Qlast]!=quadrant){
	Qlast--;
}	
//cout<<"\n\nQlast:"<<Qlast;
Qfirst=0;
while(DataID[Qfirst]!=quadrant){
	Qfirst++;
}	
//cout<<"\nQfirst:"<<Qfirst;


if(search_time<Time[Qfirst]){
	search_time=Time[Qfirst];
	i=0;
}

if(search_time>Time[Qlast]){
	search_time=Time[Qlast-5];
	i=Qlast-7;
	//cout<<"\n search time now second last time....";
	//printf(" %lf",search_time);
}
//printf("\n\nsearch time: %lf",search_time);
//cout<<"\n quadrant is "<<quadrant;

//printf("\n####start i:%ld",i);
while(flag!=0){
while(Time[i]<=search_time){
		i++;
}

if(!(Time[i]>search_time)){
	i--;
}
if(DataID[i]==quadrant){
	flag=0;
}
else{
	i++;
}

}
//cout<<"\ni: "<<i;
//commented beacause to get proper near time value
//i--;
//printf("\nsearch time: %lf",search_time);
//printf("\nstart index : %ld",*start_index);

//printf("\n found at index: %ld",i);
//printf("\n found time is: %lf",Time[i]);

//printf("\n hkvalue at found time is: %d",HKChannelNo[i]);
//*start_index=i;
found_time=i;

if(quadrant==0){
*start_indexQ0=found_time-1;
}
if(quadrant==1){
*start_indexQ1=found_time-1;
}
if(quadrant==2){
*start_indexQ2=found_time-1;
}
if(quadrant==3){
*start_indexQ3=found_time-1;
}

//printf("\n found at index: %ld",i);
//printf("\n found time is: %lf",Time[i]);
found_time_quadrant=found_time;
//printf("\ndataID is 1 at: %ld",i);
count=0;
int def=0;
if(DataID[i]==quadrant){
	while(count!=8){
		if(DataID[i]==quadrant){
			hkval=HKChannelNo[i];
			//cout<<"\nhkval: "<<hkval;
			if(adc[hkval]==-1){
				adc[hkval]=ADCOutput[i];
				count++;
			}
		}	
		i--;
		//cout<<"\ni: "<<i;
		if(i<=0){
			/*adc[0]=-1;
			adc[1]=-1;
			adc[2]=-1;
			adc[3]=-1;
			adc[4]=-1;
			adc[5]=-1;
			adc[6]=-1;
			adc[7]=-1;*/
		//cout<<"\n output array cleared..";
			//cout<<"\n\nbreaking 1st if";
			def=1;
			break;
		}
	}
}
//printf("\n\n i test: %ld",i);
if(i<=0){
count=0;
i=found_time;
while((DataID[i]!=quadrant) && (search_time!=Time[i])){
	i++;
}
//printf("\n\n i test_new: %ld",i);
if(DataID[i]==quadrant){
	while(count!=8){
		if(DataID[i]==quadrant){
			hkval=HKChannelNo[i];
			//cout<<"\nhkval: "<<hkval;
			if(adc[hkval]==-1){
				adc[hkval]=ADCOutput[i];
				count++;
			}
		}	
		i++;
		if(i>=nrowshdr){
			def=2;
		//cout<<"\n\nbreaking 2nd if";
		break;
		}
			
	}


}
}
/*if(def!=0){
	adc[0]=0;
	adc[1]=1822;
	adc[2]=704;
	adc[3]=917;
	adc[4]=1818;
	adc[5]=1813;
	adc[6]=1267;
	adc[7]=285;
}*/	

hkparam[0]=adc[0];
hkparam[1]=adc[1]/409.6;
hkparam[2]=100*(1.67-(adc[2]/409.6));
hkparam[3]=adc[3]/409.6;
hkparam[4]=adc[4]/409.6;
hkparam[5]=adc[5]/409.6;
hkparam[6]=adc[6]/409.6;
hkparam[7]=adc[7]/409.6;

*error_count=ErrorCount[found_time_quadrant];
*boot_page_no=BootPageNo[found_time_quadrant];
*alphacounter=AlphaCount[found_time_quadrant];
*vetocounter=VetoCount[found_time_quadrant];
*modeid=ModeID[found_time_quadrant];


/*cout<<"\n QUADRANT: Q"<<quadrant;
printf("\n\nhkparam[0]: %f",hkparam[0]);
printf("\nhkparam[1]: %f",hkparam[1]);
printf("\nhkparam[2]: %f",hkparam[2]);
printf("\nhkparam[3]: %f",hkparam[3]);
printf("\nhkparam[4]: %f",hkparam[4]);
printf("\nhkparam[5]: %f",hkparam[5]);
printf("\nhkparam[6]: %f",hkparam[6]);
printf("\nhkparam[7]: %f",hkparam[7]);
printf("\nerror_count: %d",*error_count);
printf("\nboot_page_no: %d",*boot_page_no);
printf("\nalphacounter: %ld",*alphacounter);
printf("\nvetocount: %ld",*vetocounter);
printf("\nmodeid: %d",*modeid);
*/

		
}


long getNumrows(const char *infile, int hdunum)
{
	long numrows=0;
	int status=0;
	fitsfile *fptr;
	fits_open_file(&fptr,infile,READONLY, &status);
	if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

	//move hdu
	fits_movabs_hdu(fptr,hdunum,NULL,&status);
	fits_get_num_rows(fptr,&numrows,&status);
	fits_close_file(fptr,&status);

	return numrows;
}


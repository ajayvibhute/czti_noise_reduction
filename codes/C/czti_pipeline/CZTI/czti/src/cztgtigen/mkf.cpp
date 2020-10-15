#include "mkf.h"

using namespace std;

void MKF_Threshold::display(){
    LOG(INFO)<< parameter << ":" << range << "  "<< LV << " to " << UV; 
}

int readMKF_Threshold(char *file,vector<MKF_Threshold> &mkf_thr){
    ifstream fptr;
    MKF_Threshold mkfparam;
    char line[100];
    string str;
    fptr.open(file,ios::in);
    if(fptr.fail()){
        LOG(ERROR)<<"Unable to read file:"<<file;
        return (EXIT_FAILURE);
    }
    int lc=0;
    int index=0;
    vector<string> words;
    fptr.clear();
    fptr.seekg(ios::beg);
    while(fptr>>str){
        //cout<<endl<<str;
        index=str.find('=',0);
        mkfparam.parameter=str.substr(0,index);
        mkfparam.range=str.substr(index+1);
        index=mkfparam.range.find('-',0);
        mkfparam.LV=atoi((mkfparam.range.substr(0,index)).c_str());
        mkfparam.UV=atoi((mkfparam.range.substr(index+1)).c_str());
        mkf_thr.push_back(mkfparam);
        lc++;
    }
    fptr.close();
    return (EXIT_SUCCESS);
    //cout<<endl<<"Line count:"<<lc<<endl;
}

int readMKF(char *file,char *mkfextname,vector<MKFparam> &mkfparam){
    fitsfile *fptr;
    int status=0;
    fits_open_file(&fptr,file,READONLY,&status);
    if(status!=0) { fits_report_error(stderr,status);  return status; }
    fits_movnam_hdu(fptr,BINARY_TBL,mkfextname,0,&status);
    if(status)  {fits_report_error(stderr,status); return (EXIT_FAILURE);  }
    long nrows=0;
    fits_get_num_rows(fptr,&nrows,&status);
    if(nrows==0){
        LOG(ERROR)<<"No data found in MKF";
        return (EXIT_FAILURE);
    }
    MKFparam param;
    for(long i=1;i<=nrows;i++){
        param.read(fptr,i);
        mkfparam.push_back(param);
    }
    fits_close_file(fptr,&status);
    if(status!=0) return status;
    return 0;
}

int MKFparam::read(fitsfile *fptr,int row){
     int status=0; 
     fits_read_col(fptr,TDOUBLE,1,row,1,1,NULL,&time,NULL,&status);  
     fits_read_col(fptr,TBYTE,2,row,1,1,NULL,&q1_modeid,NULL,&status);
     fits_read_col(fptr,TBYTE,3,row,1,1,NULL,&q1_5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,4,row,1,1,NULL,&q1_msb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,5,row,1,1,NULL,&q1_temperature,NULL,&status);
     fits_read_col(fptr,TBYTE,6,row,1,1,NULL,&q1_lsb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,7,row,1,1,NULL,&q1_2p5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,8,row,1,1,NULL,&q1_vetoHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,9,row,1,1,NULL,&q1_cztHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,10,row,1,1,NULL,&q1_msb_vetocounter,NULL,&status);
     fits_read_col(fptr,TBYTE,11,row,1,1,NULL,&q1_veto_lld,NULL,&status);
     fits_read_col(fptr,TBYTE,12,row,1,1,NULL,&q1_dvdd,NULL,&status);
     fits_read_col(fptr,TBYTE,13,row,1,1,NULL,&q1_msb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,14,row,1,1,NULL,&q1_lsb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,15,row,1,1,NULL,&q1_io9_15_8,NULL,&status);
     fits_read_col(fptr,TBYTE,16,row,1,1,NULL,&q1_io9_7_0,NULL,&status);
     fits_read_col(fptr,TBYTE,17,row,1,1,NULL,&q1_lsb_vetocounter,NULL,&status);
     
     fits_read_col(fptr,TBYTE,18,row,1,1,NULL,&q2_modeid,NULL,&status);
     fits_read_col(fptr,TBYTE,19,row,1,1,NULL,&q2_5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,20,row,1,1,NULL,&q2_msb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,21,row,1,1,NULL,&q2_temperature,NULL,&status);
     fits_read_col(fptr,TBYTE,22,row,1,1,NULL,&q2_lsb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,23,row,1,1,NULL,&q2_2p5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,24,row,1,1,NULL,&q2_vetoHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,25,row,1,1,NULL,&q2_cztHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,26,row,1,1,NULL,&q2_msb_vetocounter,NULL,&status);
     fits_read_col(fptr,TBYTE,27,row,1,1,NULL,&q2_veto_lld,NULL,&status);
     fits_read_col(fptr,TBYTE,28,row,1,1,NULL,&q2_dvdd,NULL,&status);
     fits_read_col(fptr,TBYTE,29,row,1,1,NULL,&q2_msb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,30,row,1,1,NULL,&q2_lsb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,31,row,1,1,NULL,&q2_io9_15_8,NULL,&status);
     fits_read_col(fptr,TBYTE,32,row,1,1,NULL,&q2_io9_7_0,NULL,&status);
     fits_read_col(fptr,TBYTE,33,row,1,1,NULL,&q2_lsb_vetocounter,NULL,&status);
     
     fits_read_col(fptr,TBYTE,34,row,1,1,NULL,&q3_modeid,NULL,&status);
     fits_read_col(fptr,TBYTE,35,row,1,1,NULL,&q3_5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,36,row,1,1,NULL,&q3_msb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,37,row,1,1,NULL,&q3_temperature,NULL,&status);
     fits_read_col(fptr,TBYTE,38,row,1,1,NULL,&q3_lsb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,39,row,1,1,NULL,&q3_2p5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,40,row,1,1,NULL,&q3_vetoHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,41,row,1,1,NULL,&q3_cztHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,42,row,1,1,NULL,&q3_msb_vetocounter,NULL,&status);
     fits_read_col(fptr,TBYTE,43,row,1,1,NULL,&q3_veto_lld,NULL,&status);
     fits_read_col(fptr,TBYTE,44,row,1,1,NULL,&q3_dvdd,NULL,&status);
     fits_read_col(fptr,TBYTE,45,row,1,1,NULL,&q3_msb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,46,row,1,1,NULL,&q3_lsb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,47,row,1,1,NULL,&q3_io9_15_8,NULL,&status);
     fits_read_col(fptr,TBYTE,48,row,1,1,NULL,&q3_io9_7_0,NULL,&status);
     fits_read_col(fptr,TBYTE,49,row,1,1,NULL,&q3_lsb_vetocounter,NULL,&status);
     
     fits_read_col(fptr,TBYTE,50,row,1,1,NULL,&q4_modeid,NULL,&status);
     fits_read_col(fptr,TBYTE,51,row,1,1,NULL,&q4_5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,52,row,1,1,NULL,&q4_msb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,53,row,1,1,NULL,&q4_temperature,NULL,&status);
     fits_read_col(fptr,TBYTE,54,row,1,1,NULL,&q4_lsb_czt_counter,NULL,&status);
     fits_read_col(fptr,TBYTE,55,row,1,1,NULL,&q4_2p5v_monitor,NULL,&status);
     fits_read_col(fptr,TBYTE,56,row,1,1,NULL,&q4_vetoHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,57,row,1,1,NULL,&q4_cztHVmonitor,NULL,&status);
     fits_read_col(fptr,TBYTE,58,row,1,1,NULL,&q4_msb_vetocounter,NULL,&status);
     fits_read_col(fptr,TBYTE,59,row,1,1,NULL,&q4_veto_lld,NULL,&status);
     fits_read_col(fptr,TBYTE,60,row,1,1,NULL,&q4_dvdd,NULL,&status);
     fits_read_col(fptr,TBYTE,61,row,1,1,NULL,&q4_msb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,62,row,1,1,NULL,&q4_lsb_alphacounter,NULL,&status);
     fits_read_col(fptr,TBYTE,63,row,1,1,NULL,&q4_io9_15_8,NULL,&status);
     fits_read_col(fptr,TBYTE,64,row,1,1,NULL,&q4_io9_7_0,NULL,&status);
     fits_read_col(fptr,TBYTE,65,row,1,1,NULL,&q4_lsb_vetocounter,NULL,&status);  
     
     fits_read_col(fptr,TBYTE,66,row,1,1,NULL,&syncword,NULL,&status);  
     if(status!=0) return status;
     return 0;
}
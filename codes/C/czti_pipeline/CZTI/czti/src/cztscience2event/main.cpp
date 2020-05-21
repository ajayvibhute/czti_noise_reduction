/* 
 * File: main.cpp
 * Author: Preeti, Tanul
 *
 * Created on July 6, 2012, 5:23 PM
 */

#include <cstdlib>
#include <iostream>
#include "cztscience2event.h"


using namespace std;

int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztscience2event.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);

    time_t st,et;
        checkPFILESenv();
    checkParFile("cztscience2event");
    cztscience2event par;
    par.read(argc,argv);
    par.display();
st=time(NULL);
    LOG(INFO)<<"CZTSCIENCE2EVENT STARTED....................";
            LOG(INFO) << setfill(' ') <<  setw(20) << "as";
    int status=0;
    status=par.cztscience2eventProcess();
    if(status){
        LOG(ERROR)<<"***Error in cztscience2event process***";

    }else
        LOG(INFO)<<"CZTSCIENCE2EVENT COMPLETED SUCCESSFULLY";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed---"<<et-st<<" seconds"<<endl;

    return 0;
}


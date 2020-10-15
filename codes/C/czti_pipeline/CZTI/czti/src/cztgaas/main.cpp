/* 
 * File:   main.cpp
 * Author: preeti
 *
 * Created on September 4, 2012, 3:50 PM
 */

#include <cstdlib>
#include <iostream>
#include "cztgaas_v2.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztgaas.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		

    time_t st, et;
    
    checkPFILESenv();
    checkParFile("cztgaas");

    Cztgaas param;
    param.read(argc, argv);
    param.display();
    int flag = 0;
    st = time(NULL);
	 LOG(INFO)<<"CZTGAAS STARTED....................";
    flag = param.cztgaas_process();
    if (flag) {
        LOG(ERROR) << "***Error in CZTGAAS***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << endl << "-----------CZTGAAS COMPLETED SUCCESSFULLY----------";
    et = time(NULL);
    LOG(INFO) << "Time Elapsed---" << et - st << " seconds" << endl;
    return (EXIT_SUCCESS);
}

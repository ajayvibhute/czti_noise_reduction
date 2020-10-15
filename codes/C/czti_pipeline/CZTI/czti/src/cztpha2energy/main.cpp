
#include <cstdlib>
#include "cztpha2energyv2.h"
//#include "cztpha2energy_buffer.h"
using namespace std;

/*  Program to convert PHA values in event file to PI and energy values
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztpha2energy.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);

    time_t st, et;
    checkPFILESenv();
    checkParFile("cztpha2energy");
    int flag = 0;
	
    Cztpha2energy param;

   param.read(argc, argv);
   param.display();
	st = time(NULL);
    LOG(INFO) << "CZTPHA2ENERGY STARTED.............";
    flag = param.cztpha2energy_process();
    if (flag) {
        LOG(ERROR) << "***Error in CZTPHA2ENERGY***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO) << "CZTPHA2ENERGY COMPLETED SUCCESSFULLY";
    et = time(NULL);
    LOG(INFO) << "Time Elapsed---" << et - st << " seconds" << endl;
    

	//return EXIT_SUCCESS;
}

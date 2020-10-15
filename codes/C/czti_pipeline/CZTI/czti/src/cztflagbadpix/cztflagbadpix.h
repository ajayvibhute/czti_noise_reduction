#ifndef CZTFLAGBADPIX_H
#define CZTFLAGBADPIX_H

#include "glog/logging.h"
#include "utils.h"
#include "badpixCALDB.h"

using namespace std;

class cztflagbadpix{
private:
    char modulename[NAMESIZE];
    int nbadpixFiles;
    vector <string> badpixFiles;
    char outfile[PIL_LINESIZE];
    int history;
    int clobber;
    int debug;
    
public:
    cztflagbadpix();
    void display();
    int read(int argc, char **argv);
    int read(int nbadpixFiles, vector <string> badpixFiles, char *outfile, 
             int debug=NO, int history=YES, int clobber=YES);
    
    int cztflagbadpix_process();
    
    int get_history(vector<string> &vhistory);
};









#endif /* CZTFLAGBADPIX_H */

/* 
 * @file gtiHandler.h 
 * @author Tanul Gupta
 * @date Created on Sep 24, 2015, 1:03 PM
 * @brief Handles GTI file reading
 * @details 
 * @version 0.2.0
 */

#include "utils.h"
#include "level1handler.h"
#include "errorHandler.h"

using namespace std;

struct GTIrecord{
    double tstart;
    double tstop;
};
class GTIhandler{
private:
    string gtiFilename;
    vector <GTIrecord> GTI;
public:
    GTIhandler();
    int read_gti_file(string gtiFilename="",string gtiextname="GTI");
    int write_gti_file(fitsfile *fgti,char *extname);
    
    //GETTERS
    vector<GTIrecord> get_GTI(){return GTI;}
    //SETTERS
    void set_GTI(vector <GTIrecord> GTI){this->GTI=GTI;}
    void set_GTI(vector <double> tstart, vector <double> tstop);

    // New function which will clean reduntant entries in GTI (Mithun 08/12/15)
    int cleanup_GTI();
};

//compare function to sort gti
bool min_tstart_value(GTIrecord r1, GTIrecord r2);
//to find intersection of multiple GTI record vectors.
int find_intersecting_range(vector<GTIrecord> gti1, vector<GTIrecord> gti2, vector<GTIrecord> *gtiout);

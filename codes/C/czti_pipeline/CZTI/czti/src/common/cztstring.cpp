
#include "cztstring.h"
using namespace std;

int csvreader(char* inputString, int nelems, string* outputStringArray, char* delimiter){
	string Sdata = (string) inputString;
	string temp="";
	int pos1=0,pos2=0;
	for(int j=0; j<nelems; j++){
		pos2=Sdata.find(delimiter,pos1);
		temp = Sdata.substr(pos1,pos2-pos1);
		outputStringArray[j]=temp;
		pos1=pos2+1;
	}
	return 0;
}

int stringFinder(char *inputString, char* searchString, int pos, int* nOccurrence){
	string Sdata= (string)inputString;
	string searchData = (string)searchString;
	int temp=0,count=0;
	while(temp>-1){
		temp=Sdata.find(searchData, pos);
		pos=temp+1;
		count++;
	}
	*nOccurrence = count-1;
	return 0;
}

int csvStrToInt(char *inputString, char* delimiter, int* outputIntArray, int* nelems){
	int nOccurrence;
	int nelements;
	string* outputStringArray;

	stringFinder(inputString, delimiter, 0, &nOccurrence);

	nelements = nOccurrence + 1;
	outputStringArray = new string[nelements];
	csvreader(inputString, nelements, outputStringArray,  delimiter);

	for(int i=0; i<nelements; i++){
		//cout << outputStringArray[i]<<endl;
		outputIntArray[i]= atoi(outputStringArray[i].c_str());
		//cout << outputIntArray[i]<<endl;
	}
        
                  *nelems = nelements;;

	return 0;
}


vector<int> get_quadsToProcess(char* quadsToProcess) {
    int* quadArray;
    int noccurrence;
    stringFinder(quadsToProcess, ",", 0, &noccurrence);
    int nelems_quad = noccurrence + 1;
    quadArray = new int[nelems_quad];
    vector <int> vec_quads;
    string errorMsg="";

    if(csvStrToInt(quadsToProcess, ",", quadArray, &nelems_quad)){
        errorMsg = "***Error in converting quadrant string array into integer array***"; 
        throw errorMsg; 
    }
    
    vec_quads.assign(quadArray, quadArray+nelems_quad);
    
    return vec_quads; 
}

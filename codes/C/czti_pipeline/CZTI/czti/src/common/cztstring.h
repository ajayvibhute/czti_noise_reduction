/* 
 * File:   cztstring.h
 * Author: tanul
 *
 * Created on May, 7, 2013, 02:24 PM
 */

#ifndef CZTSTRING_H
#define	CZTSTRING_H

//**********************Include files***********************
#include <iostream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>

//****************Macro variable definitions****************

//**********************************************************

using namespace std;



//*******************Function declarations*******************
/**
 * Function that reads the csv string and then stores delimited values in string array.
 * @param inputString - string provided by the user which needs to be delimited.
 * @param nelems - number of occurrences of the delimiter.
 * @param delimiter - delimiter character as to be used while delimiting the string. Ex. - "," means that comma will act as the delimiter.
 * @return outputStringArray - pointer to the ouput string array where delimited values have to be stored. 
 */
//
int csvreader(char* inputString, int nelems, string* outputStringArray, char* delimiter=",");

/**
 * function that evaluates the number of occurrence of searchString in inputString
 * @param inputString - string provided by the user which needs to be analyzed.
 * @param pos - position in the string from where search will be initiated.
 * @param searchString - string that needs to be searched in the input string. Ex."," means that comma number of occurences of comma will be searched in our input string.
 * @return nOccurence - number of occurrence of the input string. 
 */
//
int stringFinder(char *inputString, char* searchString, int pos=0, int* nOccurrence=NULL); 
//----------------------------------------------------------------

/**
 * function that evaluates the number of occurrence of searchString in inputString
 * @param inputString - string provided by the user which needs to be analyzed.
 * @param delimiter - string which will be used to delimit the input string
 * @return outputIntArray - integer array where values delimited from input string will be stored
 * @return nelems - number of elements obtained after delimiting, generally used for memory allocation.
 */
//
int csvStrToInt(char *inputString, char* delimiter, int* outputIntArray, int* nelems);

/**
 * Returns vector<int> containing number of all quadrants to be processed
 * @param quadsToProcess
 * @return 
 */
vector<int> get_quadsToProcess(char* quadsToProcess);

#endif	/* CZTSTRING_H */


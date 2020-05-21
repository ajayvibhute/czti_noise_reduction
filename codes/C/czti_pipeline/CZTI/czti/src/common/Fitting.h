   /*
 * File:   Fitting.h
 * Author: preeti
 *
 * Created on November 11, 2010, 4:15 PM
 */
/* The file contains declarations for functions required for ploynomial fitting */

#ifndef FITTING_H
#define	FITTING_H

#define MAXORDER 1
#define MINORDER 1

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "errorHandler.h"
#include "utils.h"

using namespace std;

#define SIZE 3

/* fitPolynomial: takes two x and y values in arrays x and y and their size in num.
   Returns the fitted coefficients in array coef of size order + 1*/

template<class T> int fitPolynomial(T *x,T *y,int num,double *coef,int order);
template<class T> double summation(T *series,int num,int radix);
/**
 * Fits the best polynomial between MINORDER and MAXORDER
 * @param x : Data set x
 * @param y : Data set y
 * @param num : Num of data sets
 * @param coef : Array containing coefficients
 * @param order : Best Order fitted is returned
 */
template<class T> int fitBestOrder(T *x,T *y,int num,double *coef,int *order);

template<class T> double RMSE(T *x,T *y,int num,double *coef,int order);
template<class T> T **allocateMemory(long height,long width, T val);
template<class T> void freeMemory(T **array,int h,int w);
template<class T> double calcDet(T **array,int size);
template<class T> void getNextArray(T **array,int x,int y,int size,T **nextarray);
template<class T> void transpose(T **array,int size,T **tp);
template<class T> void adjoint(T **array,int size,double **adj);
template<class T> void cofactorMatrix(T **array,int size,double **cofactor);
//for type casting array from data type X to Y
template<class X,class Y> void convert(X **array1,Y **array2,int size);
template<class T> int inverse(T **array,int size,double **inv);

void reportError(int retval);


template<class T>
void printArray(T **array,int size);



/*DEFINITIONS*/
template<class T>
T **allocateMemory(long height,long width, T val)
{ T **array;
  array=new T*[height];
  if(array==NULL){
   cout<<"\n***Unable to allocate memory***\n";
   exit(EXIT_FAILURE);
  }
  long i,j;
  for(i=0;i<height;i++){
   array[i]=new T[width];
   if(array[i]==NULL){
    cout<<"\n***Unable to allocate memory***\n";
    exit(EXIT_FAILURE);
   }
  }

  for(i=0;i<height;i++){
      for(j=0;j<width;j++){
          array[i][j]=0;
      }
  }
  return array;
}

template<class T>
void freeMemory(T **array,int h,int w){
    for(int i=0;i<h;i++)  delete[] array[i];
    delete[] array;
}

template<class T>
double calcDet(T **array,int size){
    if(size<1){
        cout<<endl<<"***Array size less than 2 not allowed***\n";
        return EXIT_FAILURE;
    }

    if(size==1){
        return array[0][0];
    }
    if(size==2){
        return (array[0][0]*array[1][1]-array[0][1]*array[1][0]);
    }

    int i,j;
    int sign;
    double det=0;
    T **nextarray,val;
    for(i=0;i<size;i++){
        if(i%2==0) sign=1;
        else sign=-1;
        nextarray=allocateMemory(size-1,size-1,val);
        getNextArray(array,0,i,size,nextarray);
        det=det+(sign*array[0][i]*calcDet(nextarray,size-1));
        //cout<<endl<<"Det:"<<(sign*array[0][i]*calcDet(nextarray,size-1));
        freeMemory(nextarray,size-1,size-1);
     }

    //cout<<"Final det:"<<det;
    return det;
}

template<class T>
void getNextArray(T **array,int x,int y,int size,T **nextarray){
    
    int i,j;

    T *temp=new T[(size-1)*(size-1)];
    if(temp==NULL) return ;

    int index=0;
    for(i=0;i<size;i++){
       for(j=0;j<size;j++){
          if(i!=x && j!=y){
             temp[index]=array[i][j];
             index++;
          }
        }
    }
    index=0;
    for(i=0;i<size-1;i++){
        for(j=0;j<size-1;j++){
            nextarray[i][j]=temp[index];
            index++;
        }
    }
    delete[] temp;
    
}

template <class T>
void transpose(T **array,int size,T **tp){
    int i,j;
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            tp[i][j]=array[j][i];
        }
    }
 
}

template<class T>
void adjoint(T **array,int size,double **adj){
    double val=1;   
    double **cofactor=allocateMemory(size,size,val);
    cofactorMatrix(array,size,cofactor);
    transpose(cofactor,size,adj);
    freeMemory(cofactor,size,size);
}

template<class T>
void cofactorMatrix(T **array,int size,double **cofactor){
    int i,j;

    T **nextarray,val;
    int sign=0;
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            nextarray=allocateMemory(size-1,size-1,val);
            getNextArray(array,i,j,size,nextarray);
            if((i+j)%2==0) sign=1;
            else sign=-1;
            cofactor[i][j]=sign*calcDet(nextarray,size-1);
            freeMemory(nextarray,size-1,size-1);
        }
    }

   
  }

template<class T>
int inverse(T **array,int size,double **inv){
    double val=1;
    double **adj=allocateMemory(size,size,val);
    adjoint(array,size,adj);
    double det=calcDet(array,size);
    if(det==0){
        cout<<"\n***Determinant is 0***\n";
        return (EXIT_FAILURE);
    }
    int i,j;
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            inv[i][j]=adj[i][j]/det;
        }
    }
    freeMemory(adj,size,size);
    return (EXIT_SUCCESS);
}

template<class X,class Y>
void convert(X **array1,Y **array2,int size){
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            array2[i][j]=(Y)array1[i][j];
        }
    }
}

//template<class T>
//int fitLinearVec(vector<T> x, vector<T> y, vector <T> slopes, vector <T> offsets){
//    ErrorHandler errHandler;
//    long nrows = x.size();
//    vector <
//    int order=1;
//    if(order>nrows-1) {
//        errHandler.severity = errERROR;
//        errHandler.errorStatus = FITTING_NOT_POSSIBLE;
//        errHandler.errorMsg = "Number of sets " + itoa(nrows) + " are less than the order of fitting " + itoa(order) + ".";
//        throw errHandler;
//    }
//      
//}
template<class T>
int fitPolynomial(T *x,T *y,int num,double *coef,int order){
    //using AX=B
    //X=(inv)A*B

    if(order>num-1){
        cout<<"\n***Could not fit polynomial of order "<<order<<" as sets are only "<<num<<"***";
        return (EXIT_FAILURE);    
    }
    double val=1;
    double **A=allocateMemory<double>(order+1,order+1,val);
    double *B=new double[order+1];

    int i,j;

    //array to store sum_x,sum_x^2,sum_x^3...etc
    double *summations=new double[order*2];

    for(i=0;i<order*2;i++)
        summations[i]=summation(x,num,i+1);

    //creating matrix A
    A[0][0]=num;
    for(i=0;i<order+1;i++){
        for(j=0;j<order+1;j++){
            if(!(i==0 && j==0))
                A[i][j]=summations[i+j-1];
        }
    }

    //cout<<"\nMatrix A:\n";
    //printArray(A,order+1);
    //cout<<endl;
    //creating matrix B
    double sum=0;
    for(i=0;i<order+1;i++){
        sum=0;
        for(j=0;j<num;j++){
            double xpy=pow(x[j],i)*y[j];
            sum+=xpy;
        }

        B[i]=sum;
    }

    //calculating inverse of A
    int flag=0;
    double **invA=allocateMemory(order+1,order+1,val);
    flag=inverse(A,order+1,invA);
    if(flag){
        cout<<"\n***Error in fitting***\n";
        return (EXIT_FAILURE);
    }
    //cout<<"\nA inverse:\n";
    //printArray(invA,order+1);
    //cout<<endl;
    
    //Multiplying inv(A) and B to get coefficients
    double tempcoef;
    for(i=0;i<order+1;i++){
        tempcoef=0;
        for(j=0;j<order+1;j++){
            tempcoef+=invA[i][j]*B[j];
        }
        coef[i]=tempcoef;
    }
    freeMemory(A,order+1,order+1);
    delete[] B;
    delete[] summations;
    return (EXIT_SUCCESS);
}

template<class T>
double RMSE(T *x,T *y,int num,double *coef,int order){

    double *derived_y=new double[num];
    double rmse=0;

    for(int i=0;i<num;i++) derived_y[i]=0;

    for(int i=0;i<num;i++){
        for(int j=0;j<order+1;j++)
            derived_y[i]=derived_y[i]+coef[j]*pow(x[i],j);
	//cout<<endl<<derived_y[i]<<"  "<<y[i];
    }

    for(int i=0;i<num;i++){
        rmse=rmse+(y[i]-derived_y[i])*(y[i]-derived_y[i]);
    }

    rmse=sqrt(rmse);
    rmse=rmse/num;
    delete[] derived_y;
    return rmse;
}

template<class T>
double summation(T *series,int num,int radix){
    double sum=0;
    int i;
    for(i=0;i<num;i++){
        sum+=pow(series[i],radix);
    } 
    return sum;
}

template<class T>
void printArray(T **array,int size){
    for(int i=0;i<size;i++){
        cout<<endl;
        for(int j=0;j<size;j++){
            cout<<setprecision(2)<<setw(12)<<array[i][j];
        }
    }
    cout<<endl;
}

template<class T> int fitBestOrder(T *x,T *y,int num,double *coef,int *order){
    int flag=0;
    double rmse[MAXORDER+1];
    double coefficients[MAXORDER+1][MAXORDER+1];
    for(int i=MINORDER;i<2;i++){
        flag=fitPolynomial(x,y,num,coefficients[i],i);
        if(flag) {
            cout<<"\n***Error in fitting for order "<<i<<"***\n";
            return (EXIT_FAILURE);
        }
        rmse[i]=RMSE(x,y,num,coefficients[i],i);
    }
    double minrmse=9999999999;
    for(int i=MINORDER;i<=num;i++){
        //cout<<"\nRMSE  "<<i<<" : "<<rmse[i];
        if(rmse[i]<minrmse){
                minrmse=rmse[i];
                *order=i;
        }
    }
    
    for(int i=0;i<=*order;i++)
        coef[i]=coefficients[*order][i];
    return (EXIT_SUCCESS);
}

template <class T> int fitLinear(T *x, T *y, int numrows, double *coef){
    int status=0;
    int i,j=0; //counter variables
    status = fitPolynomial(x,y,numrows, coef, 1);
    if(status) {
        LOG(ERROR) << "Error in Linear fitting of data.";
        return status;
    }
    
    return status;
}


#endif	/* FITTING_H */


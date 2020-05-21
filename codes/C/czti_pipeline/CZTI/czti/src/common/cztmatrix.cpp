#include "cztmatrix.h"

using namespace std;


/******************************************************************************************************
 The functions described below are related to matrix arithmetic.
 *****************************************************************************************************/
int matrixproduct(real_2d_array &A, real_2d_array &B, real_2d_array &C,
        ae_int_t m, ae_int_t n, ae_int_t p) {

    if (m == 0 || n == 0 || p == 0) {
        LOG(ERROR) << "Empty Array Found";
        return (EXIT_FAILURE);
    }
    ae_int_t i, j, k;
    int index = 0;
    double val, c;
    double *outarray;
    outarray = new double[m * p];
    if (outarray == NULL) {
        LOG(ERROR) << "Out of memory error";
        return (EXIT_FAILURE);
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            val = 0;
            for (k = 0; k < n; k++) {
                val = A(i, k) * B(k, j);
                //cout<<A(i,k)<<" * "<<B(k,j)<<" = "<<val;
                c = c + val;
            }
            outarray[index++] = c;
            if (index > (m * p)) {
                LOG(ERROR) << "***Index out of bounds error***";
                LOG(ERROR) << "Index-" << index << "  Max is- " << (m * p);
                return (EXIT_FAILURE);
            }
            c = 0;
        }
    }
    index = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < p; j++) {
            C(i, j) = outarray[index++];
            if (index > (m * p)) {
                LOG(ERROR) << "***Index out of bounds error***";
                LOG(ERROR) << "Index-" << index << "  Max is- " << (m * p);
                return (EXIT_FAILURE);
            }
        }
    delete[] outarray;
    return (EXIT_SUCCESS);
}


//size of matrix A is m x n and size of matrix B is m x p
//C will be the returned matrix of size n x p

int multiply_AinvB(float *Aarr, float *Barr, int M, int N, int P, float *C) {
    //it will use the concept
    // X=InverseA x B
    // X=Inverse(TransposeA x A) x TransposeA x B
    LOG(INFO) << "Inside multiply_AinvB function";
    LOG(INFO) << "M:" << M << "  N:" << N << "  P:" << P;
    real_2d_array A, invAtxA, At, B; //At is transpose matrix of A
    ae_int_t m = M, n = N, p = P;
    A.setlength(m, n);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            A(i, j) = Aarr[i * n + j];
    //setting values for B matrix
    B.setlength(m, p);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            B(i, j) = Barr[i * p + j];
    LOG(INFO) << "Content set for A and B ";
    At.setlength(n, m);
    rmatrixtranspose((long) m, (long) n, A, 0, 0, At, 0, 0); //taking transpose of A
    LOG(INFO) << "Transpose of A completed";
    real_2d_array AtxA; //for storing At x A

    AtxA.setlength(n, n);
    matrixproduct(At, A, AtxA, n, m, n);

    //taking inverse of AtxA using SVD
    ae_int_t uneeded = 2, vtneeded = 2, additionalmemory = 2;

    real_1d_array w;
    real_2d_array u, vt;

    w.setlength(n);
    u.setlength(n, n);
    vt.setlength(n, n);

    ae_int_t info;
    matinvreport rep;

    bool retval;
    retval = rmatrixsvd(AtxA, n, n, uneeded, vtneeded, additionalmemory, w, u, vt);
    if (retval != true) {
        LOG(ERROR) << "Error in function rmatrixsvd";
        return (EXIT_FAILURE);
    }

    real_2d_array eigenmat;
    eigenmat.setlength(w.length(), w.length());
    for (int i = 0; i < eigenmat.rows(); i++) {
        for (int j = 0; j < eigenmat.cols(); j++) {
            if (i == j) eigenmat(i, j) = w(i);
            else eigenmat(i, j) = 0;
        }
    }

    real_2d_array temp;
    temp.setlength(n, n);
    invAtxA.setlength(n, n);

    //finding inverse of A
    real_2d_array vtt, ut;
    vtt.setlength(vt.cols(), vt.rows());
    ut.setlength(u.cols(), u.rows());

    rmatrixtranspose(vt.rows(), vt.cols(), vt, 0, 0, vtt, 0, 0); //taking transpose of vt
    //rmatrixinverse(vt,info,rep);
    rmatrixinverse(eigenmat, info, rep);
    rmatrixtranspose(u.rows(), u.cols(), u, 0, 0, ut, 0, 0); //taking transpose of u
    //rmatrixinverse(u,info,rep);
    //vt x eigenmatrix = temp
    if (matrixproduct(vtt, eigenmat, temp, vt.rows(), eigenmat.rows(), eigenmat.cols())) {
        LOG(ERROR) << "***Error in matrix product of vt and eigenmat***";
        return (EXIT_FAILURE);
    }
    //temp x u = inverse of A
    if (matrixproduct(temp, ut, invAtxA, temp.rows(), u.rows(), u.cols())) {
        LOG(ERROR) << "***Error in matrix product of temp and u***";
        return (EXIT_FAILURE);
    }
    real_2d_array invA;
    invA.setlength(n, m);
    if (matrixproduct(invAtxA, At, invA, invAtxA.rows(), At.rows(), At.cols())) {
        LOG(ERROR) << "***Error in matrix product of invAtxA and At***";
        return (EXIT_FAILURE);
    }
    LOG(INFO) << "Final multiplication of Ainverse and B.........";
    //multiplying this with B
    real_2d_array invAB;
    invAB.setlength(n, p);
    if (matrixproduct(invA, B, invAB, n, m, p)) {
        LOG(ERROR) << "***Error in matrix product of invA and B***";
        return (EXIT_FAILURE);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            C[i * p + j] = invAB(i, j);
    LOG(INFO) << "Mulitplication of Ainverse and B completed";
    return (EXIT_SUCCESS);
}


/*****************************************************************************************************
 Matrix arithmetic functions over.
 ******************************************************************************************************/


#ifndef ALLOCATE_ARRAY_HPP_
#define ALLOCATE_ARRAY_HPP_

namespace array
{
    /* typedef */
    typedef double*     Double1D;
    typedef double**    Double2D;
    typedef double***   Double3D;
    typedef double****  Double4D;
    typedef double***** Double5D;

    typedef int*     Int1D;
    typedef int**    Int2D;
    typedef int***   Int3D;
    typedef int****  Int4D;
    typedef int***** Int5D;

    /* 1D Array */

    template <typename T>
    T* Allocate1dArray(const int n)
    {
        T *a = new T[n];
        return a;
    }

    template <typename T>
    void Delete1dArray(T *a)
    {
        delete[] a;
    }

    /* 2D Array*/
    template <typename T>
    T** Allocate2dArray(const int n1, const int n2)
    {
        T** a = new T*[n1];
        a[0] = new T[n1*n2];
        for (int i = 0; i < n1; ++i) {
            a[i] = a[0] + i*n2;
        }
        return a;
    }

    template <typename T>
    void Delete2dArray(T **a)
    {
        delete[] a[0];
        delete[] a;
    }

    /* 3D Array */
    template <typename T>
    T*** Allocate3dArray(const int n1, const int n2, const int n3)
    {
        T ***a = new T**[n1];
        a[0] = new T*[n1*n2];
        a[0][0] = new T[n1*n2*n3];
        for (int i = 0; i < n1; ++i) {
            a[i] = a[0] + i*n2;
            for (int j = 0; j < n2; ++j) {
                a[i][j] = a[0][0] + i*n2*n3 + j*n3;
            }
        }
        return a;
    }

    template <typename T>
    void Delete3dArray(T ***a)
    {
        delete[] a[0][0];
        delete[] a[0];
        delete[] a;
    }

    /* 4D array */
    template <typename T>
    T**** Allocate4dArray(const int n1, const int n2, const int n3, const int n4)
    {
        T ****a    = new T***[n1];
        a[0]       = new T**[n1*n2];
        a[0][0]    = new T*[n1*n2*n3];
        a[0][0][0] = new T[n1*n2*n3*n4];
        for (int i = 0; i < n1; ++i) {
            a[i] = a[0] + i*n2;
            for (int j = 0; j < n2; ++j) {
                a[i][j] = a[0][0] + i*n2*n3 + j*n3;
                for (int k = 0; k < n3; ++k) {
                    a[i][j][k] = a[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;
                }
            }
        }
        return a;
    }

    template <typename T>
    void Delete4dArray(T ****a)
    {
        delete [] a[0][0][0];
        delete [] a[0][0];
        delete [] a[0];
        delete [] a;
    }

    /* 5D array */
    template <typename T>
    T***** Allocate5dArray(const int n1, const int n2, const int n3, const int n4, const int n5)
    {
        T *****a      = new T****[n1];
        a[0]          = new T***[n1*n2];
        a[0][0]       = new T**[n1*n2*n3];
        a[0][0][0]    = new T*[n1*n2*n3*n4];
        a[0][0][0][0] = new T[n1*n2*n3*n4];
        for (int i = 0; i < n1; ++i) {
            a[i] = a[0] + i*n2;
            for (int j = 0; j < n2; ++j) {
                a[i][j] = a[0][0] + i*n2*n3 + j*n3;
                for (int k = 0; k < n3; ++k) {
                    a[i][j][k] = a[0][0][0] + i*n2*n3*n4 + j*n3*n4 + k*n4;
                    for (int l = 0; l < n4; ++l) {
                        a[i][j][k][l] = a[0][0][0][0] + i*n2*n3*n4*n5 + j*n3*n4*n5 + k*n4*n5 + l*n5;
                    }
                }
            }
        }
        return a;
    }

    template<typename T>
    void Delete5dArray(T *****a)
    {
        delete [] a[0][0][0][0];
        delete [] a[0][0][0];
        delete [] a[0][0];
        delete [] a[0];
        delete [] a;
    }

} 


#endif /* ALLOCATE_ARRAY */
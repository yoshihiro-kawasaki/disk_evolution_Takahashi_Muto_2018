#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <iostream>
#include <cmath>

/* 1D array */

template<typename T>
class Array1D
{
public:

    Array1D();
    Array1D(int n1);
    Array1D(int n1, T t);
    ~Array1D();

    Array1D(const Array1D<T>& array);
    Array1D<T> &operator= (const Array1D<T> &array);

    T &operator() (const int i) { return pdata_[i]; }
    T operator() (const int i) const { return pdata_[i]; }

    void DeleteArray();

    int GetDim1() const { return n1_; }

    void Resize(int n1);
    void Resize(int n1, const T t);

    void Fill(T t);

private:

    T *pdata_;
    int n1_;
    int size_;

};

template<typename T>
Array1D<T>::Array1D()
    : pdata_(nullptr), n1_(0), size_(0)
{

}

template<typename T>
Array1D<T>::Array1D(int n1)
    : n1_(n1)
{
    size_ = n1_;
    pdata_ = new T[size_];
}


template<typename T>
Array1D<T>::Array1D(int n1, T t)
    : n1_(n1)
{
    size_ = n1_;
    pdata_ = new T[size_];
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

template<typename T>
Array1D<T>::~Array1D()
{
    DeleteArray();
}


template<typename T>
Array1D<T>::Array1D(const Array1D<T> &array)
{
    n1_ = array.n1_;
    if (array.pdata_) {
        size_ = array.size_;
        pdata_ = new T[size_];
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
}

template<typename T>
Array1D<T> &Array1D<T>::operator= (const Array1D<T> &array)
{
    if (this != &array) {
        n1_ = array.n1_;
        size_ = array.size_;
        if (pdata_ == nullptr) {
            pdata_ = new T[size_];
        } else {
            DeleteArray();
            pdata_ = new T[size_];
        }
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
    return *this;
}

template<typename T>
void Array1D<T>::DeleteArray()
{
    if (pdata_ == nullptr) {
        return;
    } else {
        delete[] pdata_;
        pdata_ = nullptr;
        return;
    }
}

template<typename T>
void Array1D<T>::Resize(int n1)
{
    n1_ = n1;
    size_ = n1_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
        return;
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
        return;
    }
}


template<typename T>
void Array1D<T>::Resize(int n1, const T t)
{
    n1_ = n1;
    size_ = n1_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
    }

    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }

    return;
}

template<typename T>
void Array1D<T>::Fill(T t)
{
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

/* 2D array */

template<typename T>
class Array2D
{

public:

    Array2D();
    Array2D(int n1, int n2);
    Array2D(int n1, int n2, T t);
    ~Array2D();

    Array2D(const Array2D<T>& array);
    Array2D<T> &operator= (const Array2D<T> &array);

    T &operator() (const int i, const int j) { return pdata_[j + n2_*i]; }
    T operator() (const int i, const int j) const { return pdata_[j + n2_*i]; }

    void DeleteArray();

    int GetDim1() const { return n1_; }
    int GetDim2() const { return n2_; }

    void Resize(const int n1, const int n2);
    void Resize(const int n1, const int n2, const T t);

    // T* Data(int i) { return &(pdata_[n1_*i]); }

    void Fill(T t);

private:

    T *pdata_;
    int n1_;
    int n2_;
    int size_;
};

template<typename T>
Array2D<T>::Array2D()
    : pdata_(nullptr), n1_(0), n2_(0)
{

}

template<typename T>
Array2D<T>::Array2D(int n1, int n2)
{
    n1_ = n1;
    n2_ = n2;
    size_ = n1_ * n2_;
    pdata_ = new T[size_];
}

template<typename T>
Array2D<T>::Array2D(int n1, int n2, T t)
{
    n1_ = n1;
    n2_ = n2;
    size_ = n1_ * n2_;
    pdata_ = new T[size_];
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

template<typename T>
Array2D<T>::~Array2D()
{
   DeleteArray();
}

template<typename T>
Array2D<T>::Array2D(const Array2D<T> &array)
{
    n1_ = array.n1_;
    n2_ = array.n2_;
    if (array.pdata_) {
        size_ = array.size_;
        pdata_ = new T[size_];
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
}

template<typename T>
Array2D<T> &Array2D<T>::operator= (const Array2D<T> &array)
{
    if (this != &array) {
        n1_ = array.n1_;
        n2_ = array.n2_;
        size_ = array.size_;
        if (pdata_ == nullptr) {
            pdata_ = new T[size_];
        } else {
            DeleteArray();
            pdata_ = new T[size_];
        }
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
    return *this;
}

template<typename T>
void Array2D<T>::DeleteArray()
{
    if (pdata_ == nullptr) {
        return;
    } else {
        delete[] pdata_;
        pdata_ = nullptr;
        return;
    }
}

template<typename T>
void Array2D<T>::Resize(const int n1, const int n2)
{
    n1_ = n1;
    n2_ = n2;
    size_ = n1_ * n2_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
        return;
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
        return;
    }
}


template<typename T>
void Array2D<T>::Resize(const int n1, const int n2, const T t)
{
    n1_ = n1;
    n2_ = n2;
    size_ = n1_ * n2_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
    }

    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }

    return;
}

template<typename T>
void Array2D<T>::Fill(T t)
{
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}


/* 3D array */

template<typename T>
class Array3D
{

public:

    Array3D();
    Array3D(int n1, int n2, int n3);
    Array3D(int n1, int n2, int n3, T t);
    ~Array3D();

    Array3D(const Array3D<T>& array);
    Array3D<T> &operator= (const Array3D<T> &array);

    // T &operator() (const int i, const int j, const int k) { return pdata_[k + n1_*(j + n2_*i)]; }
    // T operator() (const int i, const int j, const int k) const { return pdata_[k + n1_*(j + n2_*i)]; }
    T &operator() (const int i, const int j, const int k) { return pdata_[k + n3_*(j + n2_*i)]; }
    T operator() (const int i, const int j, const int k) const { return pdata_[k + n3_*(j + n2_*i)]; }

    void DeleteArray();

    int GetDim1() const { return n1_; }
    int GetDim2() const { return n2_; }
    int GetDim3() const { return n3_;}

    void Resize(const int n1, const int n2, const int n3);
    void Resize(const int n1, const int n2, const int n3, const T t);

    void Fill(T t);

private:

    T *pdata_;
    int n1_;
    int n2_;
    int n3_;
    int size_;
};

template<typename T>
Array3D<T>::Array3D()
    : pdata_(nullptr), n1_(0), n2_(0), n3_(0)
{

}

template<typename T>
Array3D<T>::Array3D(int n1, int n2, int n3)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    size_ = n1_ * n2_ * n3_;
    pdata_ = new T[size_];
}

template<typename T>
Array3D<T>::Array3D(int n1, int n2, int n3, T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    size_ = n1_ * n2_ * n3_;
    pdata_ = new T[size_];
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

template<typename T>
Array3D<T>::~Array3D()
{
   DeleteArray();
}

template<typename T>
Array3D<T>::Array3D(const Array3D<T> &array)
{
    n1_ = array.n1_;
    n2_ = array.n2_;
    n3_ = array.n3_;
    if (array.pdata_) {
        size_ = array.size_;
        pdata_ = new T[size_];
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
}

template<typename T>
Array3D<T> &Array3D<T>::operator= (const Array3D<T> &array)
{
    if (this != &array) {
        n1_ = array.n1_;
        n2_ = array.n2_;
        n3_ = array.n3_;
        size_ = array.size_;
        if (pdata_ == nullptr) {
            pdata_ = new T[size_];
        } else {
            DeleteArray();
            pdata_ = new T[size_];
        }
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
    return *this;
}

template<typename T>
void Array3D<T>::DeleteArray()
{
    if (pdata_ == nullptr) {
        return;
    } else {
        delete[] pdata_;
        pdata_ = nullptr;
        return;
    }
}

template<typename T>
void Array3D<T>::Resize(const int n1, const int n2, const int n3)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    size_ = n1_ * n2_ * n3_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
        return;
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
        return;
    }
}


template<typename T>
void Array3D<T>::Resize(const int n1, const int n2, const int n3, const T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    size_ = n1_ * n2_ * n3_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
    }

    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }

    return;
}


template<typename T>
void Array3D<T>::Fill(T t)
{
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}


/* 4D array */

template<typename T>
class Array4D
{

public:

    Array4D();
    Array4D(int n1, int n2, int n3, int n4);
    Array4D(int n1, int n2, int n3, int n4, T t);
    ~Array4D();

    Array4D(const Array4D<T>& array);
    Array4D<T> &operator= (const Array4D<T> &array);

    T &operator() (const int i, const int j, const int k, const int l) { return pdata_[l + n4_*(k + n3_*(j + n2_*i))]; }
    T operator() (const int i, const int j, const int k, const int l) const { return pdata_[l + n4_*(k + n3_*(j + n2_*i))]; }

    void DeleteArray();

    int GetDim1() const { return n1_; }
    int GetDim2() const { return n2_; }
    int GetDim3() const { return n3_; }
    int GetDim4() const { return n4_; }

    void Resize(const int n1, const int n2, const int n3, const int n4);
    void Resize(const int n1, const int n2, const int n3, const int n4, const T t);

    void Fill(T t);

private:

    T *pdata_;
    int n1_;
    int n2_;
    int n3_;
    int n4_;
    int size_;
};

template<typename T>
Array4D<T>::Array4D()
    : pdata_(nullptr), n1_(0), n2_(0), n3_(0), n4_(0)
{

}

template<typename T>
Array4D<T>::Array4D(int n1, int n2, int n3, int n4)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    size_ = n1_ * n2_ * n3_ * n4_;
    pdata_ = new T[size_];
}

template<typename T>
Array4D<T>::Array4D(int n1, int n2, int n3, int n4, T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    size_ = n1_ * n2_ * n3_ * n4_;
    pdata_ = new T[size_];
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

template<typename T>
Array4D<T>::~Array4D()
{
   DeleteArray();
}

template<typename T>
Array4D<T>::Array4D(const Array4D<T> &array)
{
    n1_ = array.n1_;
    n2_ = array.n2_;
    n3_ = array.n3_;
    n4_ = array.n4_;
    if (array.pdata_) {
        size_ = array.size_;
        pdata_ = new T[size_];
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
}

template<typename T>
Array4D<T> &Array4D<T>::operator= (const Array4D<T> &array)
{
    if (this != &array) {
        n1_ = array.n1_;
        n2_ = array.n2_;
        n3_ = array.n3_;
        n4_ = array.n4_;
        size_ = array.size_;
        if (pdata_ == nullptr) {
            pdata_ = new T[size_];
        } else {
            DeleteArray();
            pdata_ = new T[size_];
        }
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
    return *this;
}

template<typename T>
void Array4D<T>::DeleteArray()
{
    if (pdata_ == nullptr) {
        return;
    } else {
        delete[] pdata_;
        pdata_ = nullptr;
        return;
    }
}

template<typename T>
void Array4D<T>::Resize(const int n1, const int n2, const int n3, const int n4)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    size_ = n1_ * n2_ * n3_ * n4_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
        return;
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
        return;
    }
}


template<typename T>
void Array4D<T>::Resize(const int n1, const int n2, const int n3, const int n4, const T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    size_ = n1_ * n2_ * n3_ * n4_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
    }

    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }

    return;
}

template<typename T>
void Array4D<T>::Fill(T t)
{
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}


/* 5D array */

template<typename T>
class Array5D
{

public:

    Array5D();
    Array5D(int n1, int n2, int n3, int n4, int n5);
    Array5D(int n1, int n2, int n3, int n4, int n5, T t);
    ~Array5D();

    Array5D(const Array5D<T>& array);
    Array5D<T> &operator= (const Array5D<T> &array);

    T &operator() (const int i, const int j, const int k, const int l, const int m) { return pdata_[m + n5_*(l + n4_*(k + n3_*(j + n2_*i)))]; }
    T operator() (const int i, const int j, const int k, const int l, const int m) const { return pdata_[m + n5_*(l + n4_*(k + n3_*(j + n2_*i)))]; }

    void DeleteArray();

    int GetDim1() const { return n1_; }
    int GetDim2() const { return n2_; }
    int GetDim3() const { return n3_; }
    int GetDim4() const { return n4_; }
    int GetDim5() const { return n5_; }

    void Resize(const int n1, const int n2, const int n3, const int n4, const int n5);
    void Resize(const int n1, const int n2, const int n3, const int n4, const int n5, const T t);

    void Fill(T t);

private:

    T *pdata_;
    int n1_;
    int n2_;
    int n3_;
    int n4_;
    int n5_;
    int size_;
};

template<typename T>
Array5D<T>::Array5D()
    : pdata_(nullptr), n1_(0), n2_(0), n3_(0), n4_(0), n5_(0)
{

}

template<typename T>
Array5D<T>::Array5D(int n1, int n2, int n3, int n4, int n5)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    n5_ = n5;
    size_ = n1_ * n2_ * n3_ * n4_ * n5_;
    pdata_ = new T[size_];
}

template<typename T>
Array5D<T>::Array5D(int n1, int n2, int n3, int n4, int n5, T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    n5_ = n5;
    size_ = n1_ * n2_ * n3_ * n4_;
    pdata_ = new T[size_];
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}

template<typename T>
Array5D<T>::~Array5D()
{
   DeleteArray();
}

template<typename T>
Array5D<T>::Array5D(const Array5D<T> &array)
{
    n1_ = array.n1_;
    n2_ = array.n2_;
    n3_ = array.n3_;
    n4_ = array.n4_;
    n5_ = array.n5_;
    if (array.pdata_) {
        size_ = array.size_;
        pdata_ = new T[size_];
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
}

template<typename T>
Array5D<T> &Array5D<T>::operator= (const Array5D<T> &array)
{
    if (this != &array) {
        n1_ = array.n1_;
        n2_ = array.n2_;
        n3_ = array.n3_;
        n4_ = array.n4_;
        n5_ = array.n5_;
        size_ = array.size_;
        if (pdata_ == nullptr) {
            pdata_ = new T[size_];
        } else {
            DeleteArray();
            pdata_ = new T[size_];
        }
        for (int i = 0; i < size_; ++i) {
            pdata_[i] = array.pdata_[i];
        }
    }
    return *this;
}

template<typename T>
void Array5D<T>::DeleteArray()
{
    if (pdata_ == nullptr) {
        return;
    } else {
        delete[] pdata_;
        pdata_ = nullptr;
        return;
    }
}

template<typename T>
void Array5D<T>::Resize(const int n1, const int n2, const int n3, const int n4, const int n5)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    n5_ = n5;
    size_ = n1_ * n2_ * n3_ * n4_ * n5_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
        return;
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
        return;
    }
}


template<typename T>
void Array5D<T>::Resize(const int n1, const int n2, const int n3, const int n4, const int n5, const T t)
{
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n4_ = n4;
    n5_ = n5;
    size_ = n1_ * n2_ * n3_ * n4_ * n5_;
    if (pdata_ == nullptr) {
        pdata_ = new T[size_]();
    } else {
        DeleteArray();
        pdata_ = new T[size_]();
    }

    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }

    return;
}

template<typename T>
void Array5D<T>::Fill(T t)
{
    for (int i = 0; i < size_; ++i) {
        pdata_[i] = t;
    }
}




#endif /* ARRAY_HPP */
/*
 * Declaration of the vector class.
 * It includes some global constans and functions and the definitions of
 * the following vector classes:
 * cVector for any type of vectors,
 * cOffsetVector for vectors that have a shifted index,
 * dvec for double-precission vectors, and
 * Vec3D for vectors with 3 components, from which dVec3D, fVec3D and iVec3D inherit.
 * Jonas D. De Basabe (jonas@cicese.mx)
 * File:   cvector.h
 * 2023/02/03 - Reimplemented the base object to use std::vector instead of pointers
 */
//---------------------------------------------------------------------------
#ifndef CVECTORHDR
#define CVECTORHDR

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <cassert>
#include <algorithm>
#include <vector>

//---------------------------------------------------------------------------
// A small value
const double EPSILON= std::numeric_limits<float>::epsilon();
//---------------------------------------------------------------------------
// elevate any type of number to the second power
template <typename T>
inline T sqr(T a)
{
    return a*a;
}
//---------------------------------------------------------------------------
// elevate any type of number to the third power
template <typename T>
inline T cub(T a)
{
    return a*a*a;
}
//---------------------------------------------------------------------------
// Absolute value
template <typename T>
inline T Absolute(T a)
{
    return (a<0) ? -a : a;
}
//---------------------------------------------------------------------------
// Return the square root of the absolute value of any type of number
template <typename T>
inline T sqrtabs(T a)
{
    return sqrt(Absolute(a));
}
//---------------------------------------------------------------------------
// Check if the input is NAN
template <typename T>
inline bool IsNan(T value)
{
    return (value != value);
}
//---------------------------------------------------------------------------
// Check if the input is infinity
template <typename T>
inline bool IsInf(T value)
{
    return std::numeric_limits<T>::has_infinity &&
        value == std::numeric_limits<T>::infinity();
}
//---------------------------------------------------------------------------
// Check if it is a small value
template <typename T>
inline bool is_small(T value, T tolerance=2.0)
{
    return Absolute(value) <= tolerance*std::numeric_limits<T>::epsilon();
}
//---------------------------------------------------------------------------
template <class T>
class cVector
{
protected:
    std::vector<T> d;

public:
    cVector(void)
    { }
    cVector(int m, T val=0)
    { Init(m,val); }
    cVector(cVector &v)
    {
        d= v.d;
    }
    ~cVector(void)
    {
        d.clear();
    }
    bool Init(int m)
    {
        try
        {
            if( m==0 )
                d.clear();
            else
                d.resize(m);
            if( m!=d.size() )
                throw(1);
        }
        catch(...)
        {
            return false;
        }
        return true;
    }
    bool Init(int m, T val)
    {
        int i;
        try
        {
            if( m==0 )
                d.clear();
            else
                d.resize(m);
            if( m!=d.size() )
                throw(1);
        }
        catch(...)
        {
        	return false;
        }

#pragma omp parallel for
        for( i=0; i<d.size(); ++i )
            d[i]= val;
        return true;
    }
    int size(void)
    {
        return d.size();
    }
    bool empty(void)
    {
        return d.empty();
    }
    void clear(void)
    {
        d.clear();
    }
    T Max(void)
    {
        T m=d[0];
        for( T x : d )
            if( m<x )
                m= x;
        return m;
    }
    T Min(void)
    {
        T m=d[0];
        for( T x : d )
            if( m>x )
                m= x;
        return m;
    }
    T AbsMin(void)
    {
        T m=fabs(d[0]);
        for( T x : d )
            if( m>fabs(x) )
                m= fabs(x);
        return m;
    }
    void Plus(cVector<T> &v, T a)
    {
        int i;
        if( d.size()==v.size() )
        {
#pragma omp parallel for
            for( i=0; i<d.size(); ++i )
                d[i]+= a*v[i];
        }
    }
    void Sort(void)
    {
        std::sort(d.begin(), d.end());
    }
    // Compute a vector-vector dot product
    T Mult(cVector<T> &v)
    {
        int nn = ( d.size()<=v.size() ) ? d.size() : v.size();
        T res= 0;
#pragma omp parallel for  reduction(+: res)
        for( int i=0; i<nn; ++i )
            res+= d[i]*v[i];
        return res;
    }
    T Norm2(void)
    {
        T res= 0;
#pragma omp parallel for  reduction(+: res)
        for( T x : d )
            res+= sqr(x);
        return sqrt(res);
    }
    void sv(cVector<T> &v, T a)
    {
        int i;
        if( d.size()==v.size() )
        {
#pragma omp parallel for
            for( i=0; i<d.size(); ++i )
                d[i]= a*v[i];
        }
    }
    T& operator[](int i)
    {
    	return d[i];
    }
    const T& operator[](int i) const
    {
    	return d[i];
    }
    cVector &operator=(const T &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<d.size(); ++i )
            d[i]= val;
        return *this;
    }
    cVector &operator+=(cVector &v)
    {
        int i;
        if( d.size()==v.size() )
        {
#pragma omp parallel for
            for( i=0; i<d.size(); ++i )
                d[i]+= v[i];
        }
        return *this;
    }
    cVector &operator-=(cVector &v)
    {
        int i;
        if( d.size()==v.size() )
        {
#pragma omp parallel for
            for( i=0; i<d.size(); ++i )
                d[i]-= v[i];
        }
        return *this;
    }
    cVector &operator*(const T &val)
    {
#pragma omp parallel for
        for( T x : d )
            x*= val;
        return *this;
    }
    cVector &operator*=(const T &val)
    {
#pragma omp parallel for
        for( T x : d )
            x*= val;
        return *this;
    }
    cVector &operator/=(cVector &v)
    {
        int i;
        if( d.size()==v.size() )
        {
#pragma omp parallel for
            for( T x : d )
                x/= v[i];
        }
        return *this;
    }
    bool operator==(cVector &v)
    {
        for(int i=0; i<d.size(); ++i)
            if(d[i]!=v[i])
                return false;
        return true;
    }
    bool operator>=(cVector &v)
    {
        for(int i=0; i<d.size(); ++i)
            if(d[i]<v[i])
                return false;
        return true;
    }
    bool operator<=(cVector &v)
    {
        for(int i=0; i<d.size(); ++i)
            if(d[i]>v[i])
                return false;
        return true;
    }
    bool operator>(cVector &v)
    {
        for(int i=0; i<d.size(); ++i)
            if(d[i]<=v[i])
                return false;
        return true;
    }
    bool operator<(cVector &v)
    {
        for(int i=0; i<d.size(); ++i)
            if(d[i]>=v[i])
                return false;
        return true;
    }
    bool IsNan(void)
    {
        for( T x : d )
            if( isnan(x) || isinf(x) )
                return true;
        return false;
    }
    bool IsSmall(void)
    {
            for( T x : d )
            if( !is_small(x) )
                return false;
        return true;
    }
};

template <class T>
class cOffsetVec:  public cVector<T>
{
protected:
    int offset;

public:
    cOffsetVec(void): cVector<T>()
    {offset=0;}
    cOffsetVec(int m, T val): cVector<T>(m,val)
    {offset=0;}
    cOffsetVec(int m, int Offset, T val): cVector<T>(m,val)
    {SetOffset(Offset);}
    void Init(int m, int Offset, T val)
    {cVector<T>::Init(m,val); SetOffset(Offset);}
    void SetOffset(int offs)
    { offset= offs; }
    int GetOffset(void)
    { return offset; }
    T& operator[](int i)
    {
        if(i<offset || i>=offset+cVector<T>::d.size()) exit(-1);
        return cVector<T>::d[i-offset];
    }
    const T& operator[](int i) const
    {
        if(i<offset || i>=offset+cVector<T>::d.size()) exit(-1);
        return cVector<T>::d[i-offset];
    }
    cOffsetVec &operator=(const double &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<cVector<T>::d.size(); ++i )
            cVector<T>::d[i]= val;
        return *this;
    }
};

class dvec : public cVector<double>
{
public:
    dvec(int m=0, double val=0.0) : cVector<double>(m,val)
    {}
    dvec(cVector<double> &v) : cVector<double>(v)
    {}
    dvec &operator=(const double &val)
    {
        int i;
#pragma omp parallel for
        for( double x : d )
            x= val;
        return *this;
    }
};

template <class T>
class Vec3D
{
public:
    T x,y,z;

    Vec3D(void)
    { Init(0); }

    Vec3D(T val)
    { Init(val); }

    Vec3D(T xval, T yval, T zval)
    { Init(xval, yval, zval); }

    void Init(T val)
    { x = y = z = val; }

    void Init(T xval, T yval, T zval)
    { x = xval; y = yval; z = zval; }

    void Display(void)
    {
        using namespace std;
        cout<<"("<<x<<", "<<y<<", "<<z<<")\n";
    }

    void Cross(Vec3D<T> &in1, Vec3D<T> &in2)
    {
        x= in1.y*in2.z - in1.z*in2.y;
        y= in1.z*in2.x - in1.x*in2.z;
        z= in1.x*in2.y - in1.y*in2.x;
    }

    T& operator[](int i)
    {
        switch( i )
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
        }
        return z;
    }
    const T& operator[](int i) const
    {
        switch( i )
        {
            case 0: return x;
            case 1: return y;
            case 2: return z;
        }
        //return NULL;
    }

    T Dot(Vec3D<T> &in1)
    {
        return in1.x*x + in1.y*y + in1.z*z;
    }

    T Norm(void)
    {
        return sqrt(sqr(x) + sqr(y) + sqr(z));
    }

    T &operator*(Vec3D<T> &in)
    {
        return Dot(in);
    }

    Vec3D<T> &operator+=(Vec3D<T> &in)
    {
        x+= in.x; y+= in.y; z+= in.z;
        return *this;
    }

    Vec3D<T> &operator-=(Vec3D<T> &in)
    {
        x-= in.x; y-= in.y; z-= in.z;
        return *this;
    }

    Vec3D<T> &operator*=(T val)
    {
        x*= val; y*= val; z*= val;
        return *this;
    }

    Vec3D<T> &operator/=(T val)
    {
        x/= val; y/= val; z/= val;
        return *this;
    }

    Vec3D<T> &operator=(T in)
    {
        x= in; y= in; z= in;
        return *this;
    }

    Vec3D<T> &operator=(Vec3D<T> &in)
    {
        x= in.x; y= in.y; z= in.z;
        return *this;
    }

    bool operator==(Vec3D<T> &in)
    {
        return x==in.x && y==in.y && z==in.z;
    }
};

class fVec3D : public Vec3D<float>
{
public:
    float Dist(fVec3D &in)
    {
        return sqrt(sqr(x - in.x) + sqr(y - in.y) + sqr(z - in.z));
    }
};

class dVec3D : public Vec3D<double>
{
public:
    double Dist(dVec3D &in)
    {
        return sqrt(sqr(x - in.x) + sqr(y - in.y) + sqr(z - in.z));
    }
    dVec3D &operator=(fVec3D &in)
    {
        x= in.x; y= in.y; z= in.z;
        return *this;
    }
};

class iVec3D : public Vec3D<int>
{
};

#endif

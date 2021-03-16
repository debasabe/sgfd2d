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
    bool fowndata;
    int n;
    T *d;

public:
    cVector(void)
    { n=0; d=NULL; fowndata= true; }
    cVector(int m, T val=0)
    { d=NULL; Init(m,val); }
    cVector(cVector &v)
    {
        fowndata= false; // do not free memory on destroy
        n= v.n; d= v.d;
    }
    ~cVector(void)
    { Empty(); }
    bool Init(int m)
    {
        try
        {
            fowndata= true;
            if( m==0 )
                n=0;
            else if( d==NULL )
            {
                n=m;
                d= new T[n];
                if( !d )
                	throw(1);
            }
            else if( n!=m )
            {
                delete[]d;
                n=m;
                d= new T[n];
                if( !d )
                	throw(2);
            }
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
            fowndata= true;
            if( m==0 )
                n=0;
            else if( d==NULL )
            {
                n=m;
                d= new T[n];
                if( !d )
                	throw(1);
            }
            else if( n!=m )
            {
                delete[]d;
                n=m;
                d= new T[n];
                if( !d )
                	throw(2);
            }
        }
        catch(...)
        {
        	return false;
        }

#pragma omp parallel for
        for( i=0; i<n; ++i )
            d[i]= val;
        return true;
    }
    int Dim(void)
    { return n; }
    T* Data(void)
    { return d; }
    bool IsEmpty(void)
    { return (d==NULL); }
    void Empty(void)
    { n= 0; if(fowndata) delete[]d; d=NULL; }
    void Swap(cVector &v)
    {
        int nn= v.n;
        T *daux= v.d;
        v.d= d; d= daux;
        v.n= n; n= nn;
    }
    void Copy(cVector<T> &v)
    {
        int i;
        if( n!=v.n )
            Init(v.n,0);
#pragma omp parallel for
        for( i=0; i<n; ++i )
            d[i]= v[i];
    }
    void Display(void)
    {
        for( int i=0; i<n; ++i )
            std::cout<< d[i] << " ";
        std::cout<<std::endl;
    }
    T Max(void)
    {
        T m=d[0];
        for( int i=1; i<n; ++i )
            if( m<d[i] )
                m= d[i];
        return m;
    }
    T Min(void)
    {
        T m=d[0];
        for( int i=1; i<n; ++i )
            if( m>d[i] )
                m= d[i];
        return m;
    }
    T AbsMin(void)
    {
        T m=fabs(d[0]);
        for( int i=1; i<n; ++i )
            if( m>fabs(d[i]) )
                m= fabs(d[i]);
        return m;
    }
    void Plus(cVector<T> &v, T a)
    {
        int i;
        if( n==v.n )
        {
#pragma omp parallel for
            for( i=0; i<n; ++i )
                d[i]+= a*v[i];
        }
    }
    void Sort(void)
    {
    	T aux;
        for(int i=0,j; i<n-1; ++i)
        {
            for( j=i+1; j<n; ++j )
                if(d[i]>d[j])
                {
                    aux= d[i];
                    d[i]= d[j];
                    d[j]= aux;
                }
        }
    }
    // Compute a vector-vector dot product
    T Mult(cVector<T> &v)
    {
        int nn = ( n<=v.Dim() ) ? n : v.Dim();
        T res= 0;
#pragma omp parallel for  reduction(+: res)
        for( int i=0; i<nn; ++i )
            res+= d[i]*v[i];
        return res;
    }//*/
    T Norm2(void)
    {
        T res= 0;
#pragma omp parallel for  reduction(+: res)
        for( int i=0; i<n; ++i )
            res+= sqr(d[i]);
        return sqrt(res);
    }
    void sv(cVector<T> &v, T a)
    {
        int i;
        if( n==v.n )
        {
#pragma omp parallel for
            for( i=0; i<n; ++i )
                d[i]= a*v[i];
        }
    }
    T& operator[](int i)
    {
    	assert(i>=0 && i<n );
    	return d[i];
    }
    const T& operator[](int i) const
    {
    	assert(i>=0 && i<n );
    	return d[i];
    }
    cVector &operator=(const T &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<n; ++i )
            d[i]= val;
        return *this;
    }
    cVector &operator=(cVector &v)
    {
        Empty();
        fowndata= false; // do not free memory on destroy
        n= v.n; d= v.d;
        return *this;
    }
    cVector &operator+=(cVector &v)
    {
        int i;
        if( n==v.n )
        {
#pragma omp parallel for
            for( i=0; i<n; ++i )
                d[i]+= v[i];
        }
        return *this;
    }
    cVector &operator-=(cVector &v)
    {
        int i;
        if( n==v.n )
        {
#pragma omp parallel for
            for( i=0; i<n; ++i )
                d[i]-= v[i];
        }
        return *this;
    }
    cVector &operator*(const T &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<n; ++i )
            d[i]*= val;
        return *this;
    }
    cVector &operator*=(const T &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<n; ++i )
            d[i]*= val;
        return *this;
    }
    cVector &operator/=(cVector &v)
    {
        int i;
        if( n==v.n )
        {
#pragma omp parallel for
            for( i=0; i<n; ++i )
                d[i]/= v[i];
        }
        return *this;
    }
    bool operator==(cVector &v)
    {
        for(int i=0; i<n; ++i)
            if(d[i]!=v[i])
                return false;
        return true;
    }
    bool operator>=(cVector &v)
    {
        for(int i=0; i<n; ++i)
            if(d[i]<v[i])
                return false;
        return true;
    }
    bool operator<=(cVector &v)
    {
        for(int i=0; i<n; ++i)
            if(d[i]>v[i])
                return false;
        return true;
    }
    bool operator>(cVector &v)
    {
        for(int i=0; i<n; ++i)
            if(d[i]<=v[i])
                return false;
        return true;
    }
    bool operator<(cVector &v)
    {
        for(int i=0; i<n; ++i)
            if(d[i]>=v[i])
                return false;
        return true;
    }
    bool IsNan(void)
    {
        for(int i=0; i<n; ++i)
            if( isnan(d[i]) || isinf(d[i]) )
                return true;
        return false;
    }
    bool IsSmall(void)
    {
        for(int i=0; i<n; ++i)
            if( !is_small(d[i]) )
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
    cOffsetVec(cOffsetVec &v)
    {
        cVector<T>::fowndata= false; // do not free memory on destroy
        cVector<T>::n= v.n; cVector<T>::d= v.d; offset= v.offset;
    }
    void Init(int m, int Offset, T val)
    {cVector<T>::Init(m,val); SetOffset(Offset);}
    void SetOffset(int offs)
    { offset= offs; }
    int GetOffset(void)
    { return offset; }
    T& operator[](int i)
    {
        if(i<offset || i>=offset+cVector<T>::n) exit(-1);
        return cVector<T>::d[i-offset];
    }
    const T& operator[](int i) const
    {
        if(i<offset || i>=offset+cVector<T>::n) exit(-1);
        return cVector<T>::d[i-offset];
    }
    cOffsetVec &operator=(const double &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<cVector<T>::n; ++i )
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
    dvec(int m, double *v)
    {
    	d= v;
    	n= m;
    	fowndata= false;
    }
    dvec &operator=(cVector<double> &v)
    {
        Empty();
        fowndata= false;
        cVector<double>::n= v.Dim();
        cVector<double>::d= v.Data();
        return *this;
    }
    dvec &operator=(const double &val)
    {
        int i;
#pragma omp parallel for
        for( i=0; i<cVector<double>::n; ++i )
            cVector<double>::d[i]= val;
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

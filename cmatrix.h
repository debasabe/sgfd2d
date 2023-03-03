/*
 * Declaration of the classes to manage matrices.
 * The class cMat2 uses templates and includes all the general methods.
 * The class cBandMatrix uses efficient storage for band matrices.
 * The class dmat is for double-precission matrices.
 * Jonas D. De Basabe (jonas@cicese.mx)
 * File:   cmatrix.h
 */


//---------------------------------------------------------------------------
#ifndef CMATRIXHDR
#define CMATRIXHDR
#include <fstream>
#include <iterator>
#include "cvector.h"

//---------------------------------------------------------------------------
template <class T>
class cMat2: public cVector<T>
{
protected:
	using cVector<T>::d;
	int n1, n2;
    std::vector<T*> m;

public:
	cMat2(void)
	{
		n1= n2=0;
		cVector<T>::Init(0);
	}
	cMat2(int i1, int i2)
	{
		Init(i1,i2);
	}
	cMat2(int i1, int i2, T v)
	{
		Init(i1,i2,v);
	}
	~cMat2(void)
	{
		m.clear();
	}
	bool Init(int i1, int i2)
	{
		try
		{
			if( !cVector<T>::Init(i1*i2) )
				throw(1);
			if( d.size()==0 )
				throw(1);
			n1= i1;
			n2= i2;
			m.resize(n1);
			if( m.empty() )
				throw(2);

			for(int j=0; j<n1; ++j)
				m[j]= &d[j*n2];
		}
		catch(...)
		{
			n1= n2= 0;
			m.clear();
			return false;
		}
		return true;
	}
	bool Init(int i1, int i2, T v)
	{
		try
		{
			if( !cVector<T>::Init(i1*i2,v) )
				throw(1);
			if( d.size()==0 )
				throw(1);
			n1= i1;
			n2= i2;
			m.resize(n1);
			if( m.empty() )
				throw(2);

			for(int j=0; j<n1; ++j)
				m[j]= &d[j*n2];
		}
		catch(...)
		{
			n1= n2= 0;
			m.clear();
			return false;
		}
		return true;
	}
	int Dim1(void)
	{
		return n1;
	}
	int Dim2(void)
	{
		return n2;
	}
	void clear(void)
	{
		cVector<T>::clear();
		m.clear();
		n1= n2= 0;
	}
    T*& operator[](int i)
    {
    	assert(i>=0 && i<n1 );
    	return m[i];
    }
    const T*& operator[](int i) const
    {
    	assert(i>=0 && i<n1 );
    	return m[i];
    }
    bool Swap(int i, int j)
    {
    	T* aux;
    	if( i<0 || j<0 || i>n1 || j>n1)
    		return false;
    	aux= m[i];
    	m[i]= m[j];
    	m[j]= aux;
    	return true;
    }
    void Copy(cMat2<T> &v)
    {
        int i,j;
        if( n1!=v.n1 || n2!=v.n2 )
            Init(v.n1,v.n2);

        for( i=0; i<n1; ++i )
            for( j=0; j<n2; ++j )
                m[i][j]= v[i][j];
    }
    cMat2<T>& operator=(const T &val)
    {
        for( int i=0; i<n1; ++i )
        	for( int j=0; j<n2; ++j )
        		m[i][j]= val;
        return *this;
    }
    cMat2<T>& operator+=(cMat2<T> &v)
    {
        if( n1==v.n1 && n2==v.n2 )
        {
            for( int i=0; i<n1; ++i )
                for( int j=0; j<n2; ++j )
                    m[i][j]+= v.m[i][j];
        }
        return *this;
    }
};

class dmat : public cMat2<double>
{
public:
    dmat(void): cMat2<double>()
    {}
    dmat(int m1, int m2, double val) : cMat2<double>(m1, m2, val)
    {}
    ~dmat(void) {}
    bool Init(int m1, int m2, double val)
    { return cMat2<double>::Init(m1, m2, val); }
    dmat& operator=(const double &val)
    {
        for( int i=0; i<n1; ++i )
        	for( int j=0; j<n2; ++j )
        		m[i][j]= val;
        return *this;
    }
    bool SaveTxt(char *fname)
    {
        FILE * hFile= fopen(fname,"wt");
        if( hFile!=NULL )
        {
            for( int i=0; i<n1; ++i )
            {
                for( int j=0; j<n2; ++j )
                {
                    fprintf(hFile,"%14.4e ", m[i][j]);
                }
                fprintf(hFile,"\n");
            }
            fclose(hFile);
            return true;
        }
        else return false;
    }
    bool Save(char *fname)
    {
        using namespace std;
        const int n=sizeof(float);
        std::vector<float> f(n2);
        // Don't create the file if there is no data to save
        if( empty() || n1==0 || n2==0 )
            return false;

        // Create a binary file
        ofstream hFile(fname,ios::out|ios::binary);
        if( hFile.is_open() )
        {
            for( int i=0; i<n1; ++i )
            {
#pragma omp parallel for
                for( int j=0; j<n2; ++j )
                    f[j]= (float) m[i][j];

                hFile.write((char*)&f[0],f.size()*n);
            }
            hFile.close();
            return true;
        }
        // return false in case of errors
        return false;
    }
	void Display(void)
	{
        for( int i=0; i<n1; ++i )
        {
            for( int j=0; j<n2; ++j )
				printf("%8.3f ",m[i][j]);
			printf("\n");
        }
		printf("\n");
	}
	void Display(int i1, int i2, int j1, int j2)
	{
        for( int i=i1; i<=i2; ++i )
        {
            for( int j=j1; j<=j2; ++j )
				printf("%8.3g ",m[i][j]);
			printf("\n");
        }
		printf("\n");
	}
};

//---------------------------------------------------------------------------
template <class T>
class cMat3: public cVector<T>
{
protected:
	using cVector<T>::n;
	using cVector<T>::d;
	bool fowndata;
	int n1, n2, n3;
	T*** m;

public:
	cMat3(void)
	{
		n1= n2= n3= 0;
		m= NULL;
		cVector<T>::Init(0);
	}
	cMat3(int i1, int i2, int i3)
	{
		Init(i1,i2,i3);
	}
	cMat3(int i1, int i2, int i3, T v)
	{
		Init(i1,i2, i3,v);
	}
	~cMat3(void)
	{
		Empty();
	}
	bool Init(int i1, int i2, int i3)
	{
		try
		{
			if( !cVector<T>::Init(i1*i2*i3) )
				throw(1);
			if( n==0 )
				throw(1);
			n1= i1;
			n2= i2;
			n3= i3;
			m = new T**[n1];
			if( !m )
				throw(2);

			for(int j=0; j<n1; ++j)
			{
				m[j]= new T*[n2];
				for( int k=0; k<n2; ++k )
					m[j][k]= &d[j*n2*n3 + k*n3];
			}
		}
		catch(...)
		{
			n1= n2= n3= 0;
			m= NULL;
			return false;
		}
		return true;
	}
	bool Init(int i1, int i2, int i3, T v)
	{
		try
		{
			if( !cVector<T>::Init(i1*i2*i3,v) )
				throw(1);
			if( n==0 )
				throw(1);
			n1= i1;
			n2= i2;
			n3= i3;
			m = new T**[n1];
			if( !m )
				throw(2);

			for(int j=0; j<n1; ++j)
			{
				m[j]= new T*[n2];
				for( int k=0; k<n2; ++k )
					m[j][k]= &d[j*n2*n3 + k*n3];
			}
		}
		catch(...)
		{
			n1= n2= n3= 0;
			m= NULL;
			return false;
		}
		return true;
	}
	int Dim1(void)
	{
		return n1;
	}
	int Dim2(void)
	{
		return n2;
	}
	int Dim3(void)
	{
		return n3;
	}
    void SetOwnData(bool own)
    {
        fowndata= own;
    }
	void Empty(void)
	{
		cVector<T>::Empty();
		if( fowndata )
		{
			for( int i=0; i<n1; ++i )
				delete[]m[i];
			delete[]m;
		}
		m= NULL;
		n1= n2= n3= 0;
	}
    T*& operator[](int i)
    {
    	assert(i>=0 && i<n1 );
    	return m[i];
    }
    const T*& operator[](int i) const
    {
    	assert(i>=0 && i<n1 );
    	return m[i];
    }
    bool Swap(int i, int j)
    {
    	T* aux;
    	if( i<0 || j<0 || i>n1 || j>n1)
    		return false;
    	aux= m[i];
    	m[i]= m[j];
    	m[j]= aux;
    	return true;
    }
    void Swap(cMat3<T> &v)
    {
        int naux = v.n1;
        T* daux= v.d;
        T** maux= v.m;
        v.d = d;  	d = daux;
        v.m= m;		m= maux;
        v.n1= n1; 	n1= naux;
        naux= v.n2;
        v.n2= n2; 	n2= naux;
        naux= v.n3;
        v.n3= n3; 	n3= naux;
        naux= v.n;
        v.n= n; 	n= naux;
    }
    void Copy(cMat3<T> &v)
    {
        int i,j,k;
        if( n1!=v.n1 || n2!=v.n2 || n3!=v.n3 )
            Init(v.n1,v.n2, v.n3);
        for( i=0; i<n1; ++i )
        	for( j=0; j<n2; ++j )
        		for(k=0; k<n3; ++k)
        			m[i][j][k]= v[i][j][k];
    }
    cMat3<T>& operator=(const T &val)
    {
        for( int i=0; i<n1; ++i )
        	for( int j=0; j<n2; ++j )
        		for(int k=0; k<n3; ++k)
        			m[i][j][k]= val;
        return *this;
    }
    cMat3<T>& operator=(cMat3<T> &v)
    {
        Empty();
        fowndata= false;    // do not free memory on destroy
        n= v.n; n1= v.n1; n2= v.n2; n3= v.n3;
        d= v.d;
        m= v.m;
        return *this;
    }
    cMat3<T>& operator+=(cMat3<T> &v)
    {
        if( n1==v.n1 && n2==v.n2 )
        {
            for( int i=0; i<n1; ++i )
                for( int j=0; j<n2; ++j )
            		for(int k=0; k<n3; ++k)
            			m[i][j][k]+= v.m[i][j][k];
        }
        return *this;
    }
};
#endif

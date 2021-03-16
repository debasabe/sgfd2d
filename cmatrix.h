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
#include "cvector.h"
//---------------------------------------------------------------------------
template <class T>
class cMatOld
{
protected:
    bool fowndata;
    int n1, n2;
    cVector<T> *d;

public:
    cMatOld(void)
    {
        n1=n2=0; d=NULL; fowndata= true;
    }
    cMatOld(int m1, int m2, T val)
    {
        n1=n2=0; d=NULL;
        Init(m1,m2,val);
    }
    virtual ~cMatOld(void)
    { Empty(); }
    virtual bool Init(int m1, int m2, T val)
    {
        try
        {
            fowndata= true;
            if( n1!=m1 )
            {
                Empty();
                n1=m1;
                d= new cVector<T>[n1];
                if( !d )
                	throw(1);
            }
            n2=m2;
            for( int i=0; i<n1; ++i )
                if( !d[i].Init(n2,val) )
                	throw(2);
        }
        catch(...)
        {
            printf("cMatOld ERROR: Unable to allocate a matrix of %i x %i\n\n",m1,m2);
            return false;
        }
        return true;
    }
    bool IsEmpty(void)
    { return (d==NULL); }
    int Dim1(void)
    { return n1; }
    int Dim2(void)
    { return n2; }
    void Empty(void)
    {
        if( fowndata )
            delete[]d;
        d= NULL; n1= n2= 0;
    }
    void SetOwnData(bool own)
    {
        fowndata= own;
    }
    cVector<T>& operator[](int i)
    {
    	assert(i>=0 && i<n1 );
    	return d[i];
    }
    const cVector<T>& operator[](int i) const
    {
    	assert(i>=0 && i<n1 );
    	return d[i];
    }
    void Swap(cMatOld<T> &v)
    {
        int naux = v.n1;
        cVector<T> *taux= v.d;
        v.d = d;  d = taux;
        v.n1= n1; n1= naux;
        naux= v.n2;
        v.n2= n2; n2= naux;
    }
    void Copy(cMatOld<T> &v)
    {
        int i;
        if( n1!=v.n1 || n2!=v.n2 )
            Init(v.n1,v.n2,0);
        for( i=0; i<n1; ++i )
            d[i].Copy(v[i]);
    }
    T Max(void)
    {
        T mi, m= d[0].Max();
        for( int i=1; i<n1; ++i )
        {
            mi= d[i].Max();
            if( m<mi ) m= mi;
        }
        return m;
    }
    cMatOld& operator=(const T &val)
    {
        for( int i=0; i<n1; ++i ) d[i]= val;
        return *this;
    }
    cMatOld<T>& operator=(cMatOld<T> &v)
    {
        Empty();
        fowndata= false;    // do not free memory on destroy
        n1= v.n1; n2= v.n2;
        d= v.d;
        return *this;
    }
    cMatOld<T>& operator+=(cMatOld<T> &v)
    {
        if( n1==v.n1 && n2==v.n2 )
        {
            for( int i=0; i<n1; ++i )
                for( int j=0; j<n2; ++j )
                    d[i][j]+= v.d[i][j];
        }
        return *this;
    }
};

//---------------------------------------------------------------------------
template <class T>
class cMat2: public cVector<T>
{
protected:
	using cVector<T>::n;
	using cVector<T>::d;
	using cVector<T>::fowndata;
	int n1, n2;
	T** m;

public:
	cMat2(void)
	{
		n1= n2=0;
		m= NULL;
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
		Empty();
	}
	bool Init(int i1, int i2)
	{
		try
		{
			if( !cVector<T>::Init(i1*i2) )
				throw(1);
			if( n==0 )
				throw(1);
			n1= i1;
			n2= i2;
			m = new T*[n1];
			if( !m )
				throw(2);

			for(int j=0; j<n1; ++j)
				m[j]= &d[j*n2];
		}
		catch(...)
		{
			n1= n2= 0;
			m= NULL;
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
			if( n==0 )
				throw(1);
			n1= i1;
			n2= i2;
			m = new T*[n1];
			if( !m )
				throw(2);

			for(int j=0; j<n1; ++j)
				m[j]= &d[j*n2];
		}
		catch(...)
		{
			n1= n2= 0;
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
    void SetOwnData(bool own)
    {
        fowndata= own;
    }
	void Empty(void)
	{
		cVector<T>::Empty();
		if( fowndata )
			delete[]m;
		m= NULL;
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
    void Swap(cMat2<T> &v)
    {
        int naux = v.n1;
        T* daux= v.d;
        T** maux= v.m;
        v.d = d;  	d = daux;
        v.m= m;		m= maux;
        v.n1= n1; 	n1= naux;
        naux= v.n2;
        v.n2= n2; 	n2= naux;
        naux= v.n;
        v.n= n; 	n= naux;
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
    cMat2<T>& operator=(cMat2<T> &v)
    {
        Empty();
        fowndata= false;    // do not free memory on destroy
        n= v.n; n1= v.n1; n2= v.n2;
        d= v.d;
        m= v.m;
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

//---------------------------------------------------------------------------
template <class T>
class cMat3: public cVector<T>
{
protected:
	using cVector<T>::n;
	using cVector<T>::d;
	using cVector<T>::fowndata;
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

/*---------------------------------------------------------------------------
template <class T>
class cMat3
{
protected:
    bool fowndata;
    int n1, n2, n3;
    cMat2<T> *d;

public:
    cMat3(void)
    {
        n1=n2=n3=0; d=NULL; fowndata= true;
    }
    cMat3(int m1, int m2, int m3, T val=0.0)
    {
        n1=n2=n3=0; d=NULL;
        Init(m1,m2, m3,val);
    }
    virtual ~cMat3(void)
    { Empty(); }
    virtual bool Init(int m1, int m2, int m3, T val)
    {
        try
        {
            fowndata= true;
            if( n1!=m1 )
            {
                Empty();
                n1=m1;
                d= new cMat2<T>[n1];
                if( d==NULL )
                	throw(1);
            }
            n2=m2;
            n3=m3;
            for( int i=0; i<n1; ++i )
            	if( !d[i].Init(n2, n3, val) )
            		throw(2);
        }
        catch(...)
        {
            printf("cMat3 ERROR: Unable to allocate a matrix of %i x %i x %i\n\n",m1,m2,m3);
            return false;
        }
        return true;
    }
    bool IsEmpty(void)
    { return (d==NULL); }
    int Dim1(void)
    { return n1; }
    int Dim2(void)
    { return n2; }
    int Dim3(void)
    { return n3; }
    void Empty(void)
    {
        if( fowndata )
            delete[]d;
        d= NULL; n1= n2= n3= 0;
    }
    void SetOwnData(bool own)
    {
        fowndata= own;
    }
    cMat2<T>& operator[](int i)
    { return d[i]; }
    const cMat2<T>& operator[](int i) const
    { return d[i]; }
    void Swap(cMat3<T> &v)
    {
        int naux = v.n1;
        cMat2<T> *taux= v.d;
        v.d = d;  d = taux;
        v.n1= n1; n1= naux;
        naux= v.n2;
        v.n2= n2; n2= naux;
        naux= v.n3;
        v.n3= n3; n3= naux;
    }
    void Copy(cMat3<T> &v)
    {
        int i;
        if( n1!=v.n1 || n2!=v.n2 || n3!=v.n3 )
            Init(v.n1,v.n2, v.n3,0);
        for( i=0; i<n1; ++i )
            d[i].Copy(v[i]);
    }
    T Max(void)
    {
        T mi, m= d[0].Max();
        for( int i=1; i<n1; ++i )
        {
            mi= d[i].Max();
            if( m<mi ) m= mi;
        }
        return m;
    }
    cMat3<T>& operator=(const T &val)
    {
        for( int i=0; i<n1; ++i ) d[i]= val;
        return *this;
    }
    cMat3<T>& operator=(cMat3<T> &v)
    {
        Empty();
        fowndata= false;    // do not free memory on destroy
        n1= v.n1; n2= v.n2;
        d= v.d;
        return *this;
    }
    cMat3<T>& operator+=(cMat3<T> &v)
    {
        if( n1==v.n1 && n2==v.n2 && n3==v.n3 )
        {
            for( int i=0; i<n1; ++i )
                    d[i]+= v.d[i];
        }
        return *this;
    }
};
*/

class cBandMatrix
{
protected:
    int n, band;
    cOffsetVec<double> *d;

public:
    cBandMatrix(void)
    { n=band=0; d=NULL; }
    cBandMatrix(int m, int bnd, double val)
    { d= NULL; Init(m,bnd,val); }
    virtual ~cBandMatrix(void)
    { Empty(); }
    virtual void Init(int m1, int bnd, double val)
    {
        int i,b= (bnd - 1)/2;
        if( b>= m1 -1 )
        {
            Empty();
            n=m1; band=bnd;
            for( int i=0; i<n; ++i )
                d[i].Init(n,0,val);
        }
        else
        {
            Empty();
            n=m1; band= bnd;
            d= new cOffsetVec<double>[n];
            for( i=0; i<=b; ++i )
                d[i].Init(i + 1 + b, 0, val);
            for( i=b+1; i<n - b; ++i )
                d[i].Init(bnd,i-b,val);
            for( i=n - b; i<n; ++i )
                d[i].Init(n - i + b, i-b, val);
        }
    }
    void Empty(void)
    {
        if( d!=NULL )
        {
            delete[]d; d= NULL;
        }
        n= band= 0;
    }
    void Swap(cBandMatrix &v)
    {
        cOffsetVec<double> *aux= v.d;
        v.d= d;
        d= aux;
    }
    void Mult(dvec &vi, dvec &vo)
    {
        int i, j, jmax, b= (band - 1)/2;
        vo.Init(vi.Dim(),0.0);
#pragma omp parallel for
        for( i=0; i<n; ++i )
        {
            j=(i>b)? i - b : 0;
            jmax= (i+b<n)? i+b : n;
            for( ; j<jmax; ++j )
                vo[i]+= d[i][j]*vi[j];
        }
    }
    void Mult(dvec &vi, dvec *vo)
    {
        int i, j, jmax, b= (band - 1)/2;
        vo->Init(vi.Dim(),0.0);
#pragma omp parallel for
        for( i=0; i<n; ++i )
        {
            j=(i>b)? i - b : 0;
            jmax= (i+b<n)? i+b : n;
            for( ; j<jmax; ++j )
                (*vo)[i]+= d[i][j]*vi[j];
        }
    }
    void Mult(dvec &vi, dvec &vo, double a)
    {
        int i, j, jmax, b= (band - 1)/2;
        vo.Init(vi.Dim(),0.0);
#pragma omp parallel for
        for( i=0; i<n; ++i )
        {
            j=(i>b)? i - b : 0;
            jmax= (i+b<n)? i+b : n;
            for( ; j<jmax; ++j )
                vo[i]+= a*d[i][j]*vi[j];
        }
    }
    void saxpsy(dvec &vo, double alpha, dvec &vi1, double beta, dvec &vi2)
    {
        int i, j, jmax, b= (band - 1)/2;
        vo.Init(vi2.Dim(),0.0);
#pragma omp parallel for
        for( i=0; i<n; ++i )
        {
            j=(i>b)? i - b : 0;
            jmax= (i+b<n)? i+b : n;
            for( ; j<jmax; ++j )
                vo[i]+= alpha*d[i][j]*vi1[j];
            vo[i]+= beta*vi2[i];
        }
    }
    dvec &operator*(dvec &v)
    {
        dvec *w= new dvec(v.Dim(),0.0);
        Mult(v,*w);
        return *w;
    }
    cOffsetVec<double>& operator[](int i)
    { return d[i]; }
    const cOffsetVec<double>& operator[](int i) const
    { return d[i]; }
    cBandMatrix& operator=(const double &val)
    {
        for( int i=0; i<n; ++i ) d[i]= val;
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
    void Init(dmat &v, bool own)
    {
        Empty();
        d= v.d;
        m= v.m;
        n1= v.n1;
        n2= v.n2;
        n= v.n;
        SetOwnData(own);
        v.SetOwnData(!own);
    }
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
        float f;
        // Don't create the file if there is no data to save
        if( IsEmpty() || n1==0 || n2==0 )
            return false;

        // Create a binary file
        ofstream hFile(fname,ios::out|ios::binary);
        if( hFile.is_open() )
        {
            for( int i=0; i<n1; ++i )
            {
                for( int j=0; j<n2; ++j )
                {
                    f= (float) m[i][j];
                    hFile.write((char*)&f,n);
                }
            }
            hFile.close();
            return true;
        }
        // return false in case of errors
        return false;
    }
    bool Save(char *fname, int i1, int i2, int j1, int j2)
    {
        using namespace std;
        const int n=sizeof(float);
        float f;
        // Don't create the file if there is no data to save
        if( IsEmpty() || n1==0 || n2==0 )
            return false;
        // Create a binary file
        ofstream hFile(fname,ios::out|ios::binary);
        if( hFile.is_open() )
        {
            for( int i=i1; i<=i2 && i<n1; ++i )
            {
                for( int j=j1; j<=j2 && j<n2; ++j )
                {
                    f= (float) m[i][j];
                    hFile.write((char*)&f,n);
                }
            }
            hFile.close();
            return true;
        }
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
#endif

#include <math.h>
//#include "mextools/mextools.h"
#ifndef MLABVECTOR_H
#define MLABVECTOR_H

#include <stdio.h>
//#include <iostream.h>		// Had to comment this on VC7.1 Must check if it's really needed on VC6.0

#include "mex.h"
//#ifdef _OPENMP
//#include <omp.h>
//#endif

#undef malloc
#undef realloc
#undef free
#define malloc mxMalloc
#define realloc mxRealloc
#define free mxFree
#define printf mexPrintf

// declare main class names 
class mvector;
class mmatrix;

// Simple double vector class
// The base used for indexing (zero in C, one in Fortran) is adjustable
// No automatic resizing, no bound checking, just a lightweight vector class
template<class T>
class onebased_vector {
	protected:
		T* v;
		int owner;
		const long base_index;
			
	public:
		onebased_vector(T* const V, const long Base_index = 1) : v(V-Base_index), owner(0), base_index(Base_index) {};
		onebased_vector(const long L, const long Base_index = 1) : v(0), owner(1), base_index(Base_index)
		{
			if (L > 0) {
				v = new T[L];
				v -= base_index;
			}
		}	
		~onebased_vector() { 
			if (owner) {
				v += base_index; 
				delete[] v; 
			}
		}
		
		T operator()(const long i) const { return v[i]; }
		T& operator()(const long i) { return v[i]; }
};

// Simple vector class for use in mex-files, one based indexing, so the C++ code will look very
// similar to the m-file code
//
// Can be used to access input arguments : mvector vec(prhs[0]);
// or to create output arguments : mvector out(plhs[0] = mxCreateDoubleMatrix(24, 1, mxREAL));
//
// for (long i=1; i <= 24; i++) out(i) = 0.0;

class mvector : public onebased_vector<double> {
	protected:
		const long Length;
		
		mvector(const mvector& a); 	 // give dummy copy constructor to prevent copying 
		mvector(double* V, const long length) :
			onebased_vector<double>((double*) V), Length(length) {};
				
	public:
		mvector(const mxArray* const a) : 
			onebased_vector<double>((double*) mxGetPr(a)), Length(mxGetM(a) * mxGetN(a)) {};
			
		// create a vector from a column (1..N) of a matlab M by N matrix
		mvector(const mxArray* a, const long col) : 
			onebased_vector<double>(((double*) mxGetPr(a))+ (col-1)*mxGetM(a)), Length(mxGetM(a)) {};

		mvector(const long L) : onebased_vector<double>(L), Length(L) {};
		long length() const { return Length; }
		
		typedef double* iterator;
		iterator begin() { return &(v[1]); }
		iterator end() { return &(v[Length+1]); }	// past the end (STL style)
		
		friend class mmatrix;
};

// Simple N by M double matrix with one based indexing for rows and for columns
class mmatrix : protected onebased_vector<double> {
	protected:
		const long M;
		const long N;
		
		mmatrix(const mmatrix& a);		// give dummy copy constructor to prevent copying 
		
	public:
		mmatrix(const mxArray* const a) :
			onebased_vector<double>((double*) mxGetPr(a), mxGetM(a)+1), M(mxGetM(a)), N(mxGetN(a)) {};
		mmatrix(const long m, const long n) : 
			onebased_vector<double>(m*n, m+1), M(m), N(n) {};	
			
		long getM() const { return M; }
		long getN() const { return N; }
		
		double operator()(const long i, const long j) const { return v[i+j*M]; }
		double& operator()(const long i, const long j) { return v[i+j*M]; }
};
#endif


//#include "akimaspline.h"
#ifndef AKIMA_H
#define AKIMA_H

// Class for akima interpolation
// The akima coefficients are computed when the constructor is called.
// After construction, the coefficients cannot be altered in any way

template<class InputIterator>
class akima_interpolator {
	private:
		const long n;		// number of input knots (= pairs of x and y coordinates)
		
		// knot data
		double* x;			// knot data is copied to x and y while construction
		double* y;			// so the user can savely forget about his original data vectors
						
		// spline coefficient
		double* b;
		double* c;
		double* d;

		int iflag;	// != 0 indicates error state 
		
	public:
		// be careful, X and Y must contain N elements
		akima_interpolator(const long N, InputIterator X, InputIterator Y);
		~akima_interpolator();
		
		double eval(const double u) const;		// evalute spline at u
		
		int geterr() const { return iflag; }
};

template<class InputIterator>
akima_interpolator<InputIterator>::akima_interpolator(const long N, InputIterator X, InputIterator Y)
	: n(N), x(0), y(0), b(0), c(0), d(0), iflag(0)
{
	long i;
	
	if (n < 2) {  /* no possible interpolation */
  		iflag = 1;
  		return;
	}

	// allocate memory for spline coefficients
	b = new double[n];
	c = new double[n];	
	d = new double[n];

	// allocate and copy vectors with knot data
	x = new double[n+4];	
	y = new double[n+4];

	// allocate vectors for Steigungsdata
	
	double* m = new double[n+4];	

	x += 2;		// shift pointer so that we can access x[-2], x[-1], x[0], ..., x[n], x[n+1]		
	y += 2;		// same for y
	m += 2; 	

	int ascend = 1;
	
	// copy input point data
	for (i=0; i < n; i++) {
		x[i] = *X++;
		y[i] = *Y++;
	}
	
	for (i = 1; i < n; ++i) if (x[i] <= x[i-1]) ascend = 0;
	
	if (!ascend) {
	   iflag = 2;	// data is not given with ascending x values
	   return;
	}

	// extrapolate 2 points left and two points right
	x[-1] = x[0] + x[1] - x[2];
	x[-2] = x[-1] + x[0] - x[1];

	x[n] = x[n-1] + x[n-2] - x[n-3];
	x[n+1] = x[n] + x[n-1] - x[n-2];

	for (i=0; i < n-1; i++)
		m[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);

   
	y[-1] =(x[0]-x[-1])*(m[1]-2*m[0])+y[0];
	m[-1] =(y[0]-y[-1])/(x[0]-x[-1]);
	
	y[-2] =(x[-1]-x[-2])*(m[0]-2*m[-1])+y[-1];
	m[-2] =(y[-1]-y[-2])/(x[-1]-x[-2]);

	y[n] = (2*m[n-2]-m[n-3])*(x[n]-x[n-1])+y[n-1];
	m[n-1] = (y[n]-y[n-1])/(x[n]-x[n-1]);
	
	y[n+1] = (2*m[n-1]-m[n-2])*(x[n+1]-x[n])+y[n];	 
	m[n] = (y[n+1]-y[n])/(x[n+1]-x[n]);
	
	for (i=0; i < n; i++) {
		const double term1 = fabs(m[i+1] - m[i]);
		const double term2 = fabs(m[i-1] - m[i-2]);
		const double t1_t2 = term1 + term2;
		
		if (t1_t2 > 0)
			b[i] = (term1 * m[i-1] + term2 * m[i]) / t1_t2;
		else
			b[i] = 0;
	}
	
	// now calcualte the coefficients for the n-1 cubic polynomials
	for (i=0; i < n-1; i++) {		
		double t = x[i+1]-x[i];
	  	c[i] =(3*m[i]-2*b[i]-b[i+1]) / t;
 		d[i] =(b[i]+b[i+1]-2*m[i]) / (t * t);
	}
	
	m -= 2;
	delete[] m;
} 

template<class InputIterator>
akima_interpolator<InputIterator>::~akima_interpolator(){
	x -= 2;
	y -= 2;

	delete[] x;
	delete[] y;
	
	delete[] b;
	delete[] c;
	delete[] d;
}

template<class InputIterator>
double akima_interpolator<InputIterator>::eval(const double u) const {
	static long last;
	long  i, j, k;
	double w;
	
	i = last;
	
	if (i >= n-1) i = 0;
	if (i < 0)  i = 0;

	// perform a binary search
	if ((x[i] > u) || (x[i+1] < u)) {  
		i = 0;
		j = n-1;
		do {
			k = (i + j) / 2;         // split the domain to search 
			if (u < x[k])  j = k;    // move the upper bound 
			if (u >= x[k]) i = k;    // move the lower bound 
		} while (j > i+1);           //  there are no more segments to search 
	}
	
	last = i;

	// printf("%ld\n", i);

	w = u - x[i];
	w = y[i] + w * (b[i] + w * (c[i] + w * d[i]));
	
	return w;
}
#endif



void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[]) {       

	if (nrhs < 3) {
		mexPrintf("Akima interpolation : YY = AKIMA(X,Y,XX)");
		mexPrintf("or (insitu YY) : AKIMA(X,Y,XX, YY)");
		if (nrhs >= 1)
			mexErrMsgTxt("\n");
		else
			return;
	}
	
	/* handle matrix I/O */
	
	const double *x = (double *)mxGetPr(prhs[0]);
	const long Nx = mxGetM(prhs[0]) * mxGetN(prhs[0]); 	
	
	const double *y = (double *)mxGetPr(prhs[1]);
	const long Ny = mxGetM(prhs[1]) * mxGetN(prhs[1]); 		
	
	const double *xx = (double *)mxGetPr(prhs[2]);
	const long Nxx = mxGetM(prhs[2]) * mxGetN(prhs[2]); 		

	double *out;
	
	if (Nx < 2)
		mexErrMsgTxt("There should be at least two data points.");
	
	if (Ny != Nx)
		mexErrMsgTxt("Abscissa and ordinate vector should be of the same length.");

	if (nrhs == 3) {	
		if (mxGetM(prhs[2]) == 1)	// If row vector in - row vector out
			plhs[0] = mxCreateDoubleMatrix(1, Nxx, mxREAL);      
		else
			plhs[0] = mxCreateDoubleMatrix(Nxx, 1, mxREAL);      

		out = (double *) mxGetPr(plhs[0]);   
	}
	else {		// Using vector transmitted in input (in multiple calls avoids repeated memory alloc requests)
		const long Nyy = mxGetM(prhs[3]) * mxGetN(prhs[3]); 		
		if (Nyy != Nxx)
			mexErrMsgTxt("Abscissa and ordinate vectors (XX & YY) should be of the same length.");

		out = (double *) mxGetPr(prhs[3]);   
	}

 	akima_interpolator<const double*> akima(Nx, (const double*) x, (const double*) y);
		
	if (akima.geterr() == 1)
		mexErrMsgTxt("Need more input data points");
	else if (akima.geterr() == 2)
		mexErrMsgTxt("Data point must be given with ascending x values");
	else if (akima.geterr())
		mexErrMsgTxt("Error computing Akima coefficients");

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
	for (long i = 0; i < Nxx; i++)
		out[i] = akima.eval(xx[i]);	
}

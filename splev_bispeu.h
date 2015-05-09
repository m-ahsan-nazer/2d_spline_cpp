#ifndef SPLEV_BISPEU
#define SPLEV_BISPEU

#ifdef __cplusplus
extern"C" {
#endif
void splev_cpp_(double *t,int *n,double *c,int *k,double *x,double *y,
                int *m,int *ier);

void bispeu_cpp_(double *tx,int *nx,double *ty,int *ny,double *c,int *nc,
                 int *kx,int *ky,double *x,double *y,double *z,int *m,
                 double *wrk,int *lwrk,int *ier);
#ifdef __cplusplus
}
#endif

struct spline_1d_param
{	
/**
   Used with scipy splines
**/
	int n;
	double *t;
	double *c;
	int k;
	void operator()(int n_)
	{
		t = new double[n_];
		c = new double[n_];
	}
	~spline_1d_param()
	{
		delete []t;
		delete []c;
	}
	
};
struct spline_2d_param
{	
/**
   Used with scipy splines
**/
	int nx;
	double *tx;
	int ny;
	double *ty;
	int nc;
	double *c;
	int kx;
	int ky;
	int lwrk;
	
	void operator()(int nx_, int ny_, int nc_, int kx_, int ky_)
	{
		kx = kx_;
		ky = ky_;
		lwrk = kx+ky+2;//See Subroutine subroutine bispeu_cpp
		nx = nx_;
		ny = ny_;
		nc = nc_;
		tx = new double[nx];
		ty = new double[ny];
		c  = new double[nc];
	}
	~spline_2d_param()
	{
		delete []tx;
		delete []ty;
		delete []c;
	}
	
};

#endif

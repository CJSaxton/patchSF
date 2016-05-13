/*****************************************************************************
 *  Calculate the Kolmogorov structure function (at given order) of discrete
 *  time series data fron a 4-column input file.  Write a 4-column output
 *  file containing timescales, SF denominator, SF numerator, uncertainty.
 *
 *  The first command-line argument is the input data file.
 *
 *
 *  INPUT FILE FORMAT:
 *
 *  The input file has four columns:
 *
 *      t       = time of observation (temperal bin centroid)
 *      dt      = width of the temporal bin (or time uncertainty)  
 *      z       = signal or count-rate during this time bin
 *      dz      = uncertainty in the signal or count-rate
 *
 *
 *  OUTPUT FILE FORMAT:
 *
 *  The output file has four columns:
 *
 *      tau     = timescale of SF evaluation                    
 *      W       = a weight or normalisation denominator in the SF
 *      W S     = numerator in the structure function at tau
 *      W dS    = uncertainty propagated to the product (W S)
 *
 *  Thus the order-n structure function at timescale tau, and its uncertainty
 *  are:
 *
 *      S_n(tau)  = (W S) / W   =  (column 3) / (column 2)
 *      dS_n(tau) = (W dS) / W  = (column 4) / (column 2)
 *
 *  The code writes the denominator "W" as a separate column as a precaution.
 *  With temporally patchy data, in some numerical methods, there are sometimes
 *  circumstances when W=0 at some timescale tau.  Keeping "W" distinct enables
 *  an informed user to discard those rows from the output file.
 *
 *  The sequence of tau values is the union of a linearly spaced sequence
 *  (tmin<=tau<=tmax) and a logarithmically spaced sequence in the same domain.
 *  To increase or decrease the resolution of this output grid, vary the "try"
 *  command-line parameter.
 *
 *
 *  EXAMPLE:
 *
 *  ./patchsf signal.dat start=1. end=20. n=2 tmin=10. tmax=1.e3 ^results.dat
 *
 *
 *  COMMAND-LINE PARAMETERS & OPTIONS:
 *
 *	start	= use input data after this time
 *	end	= use input data before this time
 *	-logx	flag: take logarithm of independent variable (time)
 *	-logy	flag: take logarithm of dependent variable (signal)
 *
 * 	n	= order of the structure function
 *	try	= number of tau grid points -1
 *	tmin	= shortest time-scale (tmin <= tau)
 *	tmax	= longest time-scale (tau <= tmax)
 *
 *  Alternatively, calculate a Leahy normalised power spectrum.
 *
 *	-pds	flag: compute power density spectrum, not structure function
 *
 *  Results are written into a file named after the host machine + ".dat".
 *  An alternative output file can be named as:
 *
 *	^blah.dat	send output to "blah.dat" 
 *
 *
 *  COMPILATION:
 *
 *  	gcc -lm -lgsl -lgslcblas patchsf.c -o patchsf
 *
 *  This code depends on the installation of Gnu Scientific Library (GSL):
 *  	http://www.gnu.org/software/gsl/
 *
 *
 *
 *  HISTORY:
 *
 *	Written by Curtis J. Saxton (2011,2012).
 *	Tidying by Curtis J. Saxton (20160306).
 *
 *	Cite reference:
 *
 *	"Long-term X-ray variability of Swift J1644+57"
 *	 Curtis J. Saxton,  Roberto Soria,  Kinwah Wu,  N. Paul M. Kuin, 
 *	 2012, Monthly Notices Royal Astronomical Society 422, 1625-1639.
 *	 http://arxiv.org/abs/1201.5210
 *	 http://dx.doi.org/10.1111/j.1365-2966.2012.20739.x
 *	 DOI: 10.1111/j.1365-2966.2012.20739.x
 *
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

#include "arrays.c"

#define	TRUE	(1==1)
#define	FALSE	(1==0)
const double	pi = M_PI;
#define	max(a,b)	(((a)>(b))?a:b)
#define	min(a,b)	(((a)<(b))?a:b)

//  Global variables for the random number generator.
const gsl_rng_type * T_rng;
gsl_rng * r_rng;
double	ugauss()	{ double val=gsl_ran_ugaussian(r_rng);	return(val); }
double	uflat()	{ double val=gsl_ran_flat(r_rng,0.,1.);	return(val); }

//  Return the current time, as an integer.
unsigned long int	tnow()
{
	int	th,tm,ts;	time_t	t;	struct tm	*tbits;
	(void) time(&t);	tbits=gmtime(&t);
	th=tbits->tm_hour;	tm=tbits->tm_min;	ts=tbits->tm_sec;
	return(ts+60*(tm+60*th));
}

char	hostname[128];		// name of the host computer
char	wholecmd[512];		// command line, all options
char	profname[256];		// output profile filename

int	sfn=2;			// order of the structure function
double	*tau;			// sequence of tau times for SF calculation
double	tmin=FLT_MAX;		// lower limit of the tau domain
double	tmax=-FLT_MAX;		// upper limit of the tau domain
double	*sf_bot;		// denominator of the structure function
double	*sf_top;		// numerator of the structure function
double	*sf_unc;		// uncertainty in the numerator
double	tstart=0.;		// omit data before this time
double	tend=DBL_MAX;		// omit data after this time
int	ntries=2000;		// resolution of the tau output array

long double	*dat_t;		// data array: observation time
long double	*dat_dt;	// data array: time uncertainty / bin size
long double	*dat_z;		// data array: signal value
long double	*dat_dz;	// data array: signal uncertainty
int	ndata=0;		// number of data points intput

int	flag_pds=FALSE;		// Leahy power spectrum instead of Kolmogorov SF
int	flag_logx=FALSE;	// x values are logarithmic?
int	flag_logy=FALSE;	// y values are logarithmic?
int	mute_setup=FALSE;	// mute some setup output (unused)
int	flag_rel=FALSE;		// relative values (unused)
int	flag_rnd=FALSE;		// introduce random noise (unused)
double	aguess=1.e-5;		// amplitude of noise (unused)
int	nx=0;			// global arrya dimension (unused)
int	ny=0;			// global arrya dimension (unused)

#include "parse.c"


/****************************************************************************
 *  Robust root function, at integer power "p".
 ***************************************************************************/
double	nroot(double x, int p)
{
	if(x==0.) return(0.);
	if(x<0.) return(-pow(fabs(x),1./p)); else return(pow(fabs(x),1./p));
}


/****************************************************************************
 *  Compute and store a sequence of binomial coefficients.
 ***************************************************************************/
void	binco(int n, int *b)
{
	unsigned int	k;
	b[0]=b[n]=1;
	for(k=1;k<n;k++) b[k]=gsl_sf_choose(n,k);
}


/****************************************************************************
 *  Find and return the minimum value on a 2D array.
 ***************************************************************************/
double	findmin(double **a, int nx, int ny, int **m)
{
	int	i,j;	double	val=DBL_MAX;
	for(i=0;i<nx;i++) for(j=0;j<ny;j++) if(m[i][j]!=0)
		val=(val<a[i][j])?val:a[i][j];
	return(val);
}


/****************************************************************************
 *  Find and return the maximum value on a 2D array.
 ***************************************************************************/
double	findmax(double **a, int nx, int ny, int **m)
{
	int	i,j;	double	val=DBL_MIN;
	for(i=0;i<nx;i++) for(j=0;j<ny;j++) if(m[i][j]!=0)
		val=(val>a[i][j])?val:a[i][j];
	return(val);
}


/****************************************************************************
 *  Compute the Leahy normalised power spectra from discrete data.
 ***************************************************************************/
void	zigpds(double tau, int norder, int ndata, 
	long double w_t[], long double w_dt[],
	long double w_z[], long double w_dz[],
	double	*ans1, double *ans2, double *ans3)
{
	long double	*a,*b,*t,*dt,*z,*dz;
	long double	*Cr,*Ci,*Fr,*Fi,*Ur,*Ui,*Pr,*Pi;
	long double	*e1r,*e1i,*e2r,*e2i;
	long double	p2=2.L*pi,ft_r=0,ft_i=0,ftu_r=0,ftu_i=0,leahy=0;
	int	i,j;

	//  Copy the data into working arrays.
	t=ldvector(0,ndata);		dt=ldvector(0,ndata);
	z=ldvector(0,ndata);		dz=ldvector(0,ndata);
	for(i=0;i<ndata;i++) {
		t[i]=w_t[i];		dt[i]=w_dt[i];
		//z[i]=w_z[i];		dz[i]=w_dz[i];
		t[i]=dat_t[i];		dt[i]=dat_dt[i];
		z[i]=dat_z[i];		dz[i]=dat_dz[i];
	}

	//  Optional logarithmic modes.
	if(flag_logx) {
		for(i=0;i<ndata;i++) dt[i]=dt[i]/t[i];
		for(i=0;i<ndata;i++) t[i]=log(t[i]);
	}
	if(flag_logy) {
		for(i=0;i<ndata;i++) dz[i]=dz[i]/z[i];
		for(i=0;i<ndata;i++) z[i]=log(z[i]);
	}

	//  Compute the gradients of the line segments.
	a=ldvector(0,ndata);		b=ldvector(0,ndata);
	Cr=ldvector(0,ndata);		Ci=ldvector(0,ndata);
	Fr=ldvector(0,ndata);		Fi=ldvector(0,ndata);
	Pr=ldvector(0,ndata);		Pi=ldvector(0,ndata);
	Ur=ldvector(0,ndata);		Ui=ldvector(0,ndata);
	e1r=ldvector(0,ndata);		e1i=ldvector(0,ndata);
	e2r=ldvector(0,ndata);		e2i=ldvector(0,ndata);
	for(i=0;i<ndata-1;i++) {
		a[i]=(z[i+1]-z[i])/(t[i+1]+t[i]);
		b[i]=z[i]-a[i]*t[i];
		Cr[i]= -tau*a[i]/p2;
		Ci[i]= a[i]*t[i]-b[i];
		leahy += (dz[i+1]*dz[i+1]+dz[i]*dz[i])*(t[i+1]-t[i])/2.L;
	}	
	a[ndata-1]=a[ndata-2];		b[ndata-1]=b[ndata-2];
	Cr[ndata-1]=Cr[ndata-2];	Ci[ndata-1]=Ci[ndata-2];
	Fr[ndata-1]=Fr[ndata-2];	Fi[ndata-1]=Fi[ndata-2];
	Pr[ndata-1]=Fr[ndata-2];	Pi[ndata-1]=Fi[ndata-2];
	e1r[ndata-1]=e1r[ndata-2];	e1i[ndata-1]=e1i[ndata-2];
	e2r[ndata-1]=e2r[ndata-2];	e2i[ndata-1]=e2i[ndata-2];

	//  Compute the FT of each line interval.
	for(i=0;i<ndata-1;i++) {
		e1r[i]= cosl(p2*t[i+1]/tau)-cosl(p2*t[i]/tau);
		e1i[i]= sinl(p2*t[i+1]/tau)-sinl(p2*t[i]/tau);
		e2r[i]= sinl(p2*t[i+1]/tau);
		e2i[i]= cosl(p2*t[i+1]/tau);
		Fr[i]= Cr[i]*e1r[i]-Ci[i]*e1i[i]-a[i]*(t[i+1]-t[i])*e2i[i];
		Fi[i]= Cr[i]*e1i[i]+Ci[i]*e1r[i]+a[i]*(t[i+1]-t[i])*e2r[i];
		Fr[i] *= -tau/p2;
		Fi[i] *= -tau/p2;
		Pr[i]=-tau/p2*e1r[i]-e1i[i]*t[i]-e2i[i]*(t[i+1]-t[i]);
		Pi[i]=-tau/p2*e1i[i]+e1r[i]*t[i]+e2r[i]*(t[i+1]-t[i]);
	}
	Pr[ndata-1]=Pi[ndata-1]=0.;
	Fr[ndata-1]=Fi[ndata-1]=0.;

	//  Compute uncertainties using partial derivaties wrt observables.
	for(j=0;j<ndata;j++) Ur[j]=Ui[j]=0.;
	for(j=0;j<ndata;j++) {
		if(j>0) {
			Ur[j]= (Pr[j-1]-e1i[j-1]*t[j-1])/(t[j]-t[j-1]);
			Ui[j]= (Pi[j-1]+e1r[j-1]*t[j-1])/(t[j]-t[j-1]);
		} 
		if(j<ndata-1) {
			Ur[j] += (-Pr[j]+e1i[j]*t[j+1])/(t[j+1]-t[j]);
			Ui[j] += (-Pi[j]-e1r[j]*t[j+1])/(t[j+1]-t[j]);
		}
	}
	for(j=0;j<ndata;j++) {
		Ur[j] *= -tau/p2;
		Ui[j] *= -tau/p2;
	}

	//  Set the total amplitude^2, its uncertainty, and the denominator.
	ftu_r = ftu_i = 0.;
	for(ft_r=0.,i=0;i<ndata-1;i++) ft_r += Fr[i];
	for(ft_i=0.,i=0;i<ndata-1;i++) ft_i += Fi[i];
	for(i=0;i<ndata-1;i++) ftu_r += Ur[i]*Ur[i]*dz[i]*dz[i];
	for(i=0;i<ndata-1;i++) ftu_i += Ui[i]*Ui[i]*dz[i]*dz[i];
	ftu_r = sqrt( ftu_r );
	ftu_i = sqrt( ftu_i );
	*ans1 = leahy/2.;
	*ans2 = ft_r*ft_r + ft_i*ft_i;
	*ans3 = ftu_r*ftu_r + ftu_i*ftu_i;
	*ans3 = 2.*sqrt( ft_r*ft_r*ftu_r*ftu_r + ft_i*ft_i*ftu_i*ftu_i );
	//printf("%e\t% Le\t% Le\t% le\n",tau,ft_r,ft_i,(*ans3));
	//*ans2 = (*ans2)*(*ans2);	*ans3 = (*ans3)*(*ans3);
	printf("%e\t% Le\t% Le\t% Le\n",tau,leahy,(*ans2)/leahy,(*ans3)/leahy);

	free_ldvector(Cr,0,ndata);	free_ldvector(Ci,0,ndata);
	free_ldvector(Fr,0,ndata);	free_ldvector(Fi,0,ndata);
	free_ldvector(Pr,0,ndata);	free_ldvector(Pi,0,ndata);
	free_ldvector(Ur,0,ndata);	free_ldvector(Ui,0,ndata);
	free_ldvector(e1r,0,ndata);	free_ldvector(e1i,0,ndata);
	free_ldvector(e2r,0,ndata);	free_ldvector(e2i,0,ndata);
	free_ldvector(a,0,ndata);	free_ldvector(b,0,ndata);
	free_ldvector(t,0,ndata);	free_ldvector(dt,0,ndata);
	free_ldvector(z,0,ndata);	free_ldvector(dz,0,ndata);
}


/****************************************************************************
 * Calculate the structure function at specified order, from discrete data.
 *
 * 	tau	= evaluate the structure function at this lag-time
 * 	norder	= order of the structure function
 * 	ndata	= number at observational data
 * 	w_t	= data times
 * 	w_dt	= temporal bin width or uncertainty
 * 	w_z	= signal values at given times
 * 	w_dz	= uncertainties in signal values
 * 	ans1	= (output) denominator of the structure function at tau
 * 	ans2	= (output) numerator of the structure function at tau
 * 	ans3	= (output) uncertainty in the numerator of structure function
 *
 ****************************************************************************/
void	zigzag(double tau, int norder, int ndata, 
	long double w_t[], long double w_dt[],
	long double w_z[], long double w_dz[],
	double	*ans1, double *ans2, double *ans3)
{
	long double	*dSdz,sf_unc,sf_top,sf_bot;
	double	*a,*b,*t,*dt,*z,*dz,*tlo,*thi,*mask,*mu,*nu,*mup,*nup;
	double	*big,*bigA,*bigB,*bigG,*bigH,*C1,*C2,*C3,*C4,wnow;
	int	i,j,k,l,n1,m,*idi,*idj,*bc;

	//  Abbreviate the structure function order; binomial coefficients.
	n1=norder+1;
	bc=ivector(0,norder+1);
	binco(norder,bc);
	//for(i=0;i<=norder;i++) printf("%i\t",bc[i]);	printf("\n");

	//  Copy the data into working arrays.
	t=dvector(0,ndata);		dt=dvector(0,ndata);
	z=dvector(0,ndata);		dz=dvector(0,ndata);
	for(i=0;i<ndata;i++) {
		t[i]=w_t[i];		dt[i]=w_dt[i];
		//z[i]=w_z[i];		dz[i]=w_dz[i];
		t[i]=dat_t[i];		dt[i]=dat_dt[i];
		z[i]=dat_z[i];		dz[i]=dat_dz[i];
	}

	//  Optional logarithmic computations.
	if(flag_logx) {
		for(i=0;i<ndata;i++) dt[i]=dt[i]/t[i];
		for(i=0;i<ndata;i++) t[i]=log(t[i]);
	}
	if(flag_logy) {
		for(i=0;i<ndata;i++) dz[i]=dz[i]/z[i];
		for(i=0;i<ndata;i++) z[i]=log(z[i]);
	}

	//  Compute the gradients of the line segments.
	a=dvector(0,ndata);		b=dvector(0,ndata);
	for(i=0;i<ndata-1;i++) {
		a[i]=(z[i+1]-z[i])/(t[i+1]+t[i]);
		b[i]=z[i]-a[i]*t[i];
	}	
	a[ndata-1]=a[ndata-2];		b[ndata-1]=b[ndata-2];
	//a[ndata-1]=b[ndata-1]=0.;

	//  Count the intervals that overlap others.
	for(m=i=0;i<ndata;i++) for(j=0;j<ndata;j++)
		if(gsl_min(t[j+1]+tau,t[i+1])>gsl_max(t[j]+tau,t[i])) m++;

	//  Compute the overlap time intervals.
	idi=ivector(0,m);	idj=ivector(0,m);	mask=dvector(0,m);
	tlo=dvector(0,m);	thi=dvector(0,m);
	for(k=i=0;i<ndata;i++) for(j=0;j<ndata;j++)  {
		wnow=gsl_min(t[j+1]+tau,t[i+1]) - gsl_max(t[j]+tau,t[i]);
		if(!(wnow>0.)) continue;
		thi[k]=gsl_min(t[j+1]+tau,t[i+1]);
		tlo[k]=gsl_max(t[j]+tau,t[i]);
		idi[k]=i;	idj[k]=j;	mask[k]=thi[k]-tlo[k];
		k++;
	}

	//  Compute fractional intervals of the overlapping lines.
	mu=dvector(0,m);	 mup=dvector(0,m);
	nu=dvector(0,m);	 nup=dvector(0,m);
	for(k=0;k<m;k++) { i=idi[k];	j=idj[k];
		mu[k]=nu[k]=mup[k]=nup[k]=0.;
		if(i<ndata-1) {
			mu[k]=(tlo[k]-t[i])/(t[i+1]-t[i]);
			nu[k]=(thi[k]-t[i])/(t[i+1]-t[i]);
		}
		if(j<ndata-1) {
			mup[k]=(tlo[k]-t[j]-tau)/(t[j+1]-t[j]);
			nup[k]=(thi[k]-t[j]-tau)/(t[j+1]-t[j]);
		}
	}

	/*
	printf("%e < mu  <%e\n",findmin(mu,i,i,mask0),findmax(mu,i,i,mask0));
	printf("%e < nu  <%e\n",findmin(nu,i,i,mask0),findmax(nu,i,i,mask0));
	printf("%e < mup <%e\n",findmin(mup,i,i,mask0),findmax(mup,i,i,mask0));
	printf("%e < nup <%e\n",findmin(nup,i,i,mask0),findmax(nup,i,i,mask0));
	*/
	free_dvector(tlo,0,m);			free_dvector(thi,0,m);

	/*  Terms in the SF contributions from each line segment.  */
	big=dvector(0,m);	bigA=dvector(0,m);	bigB=dvector(0,m);
	bigG=dvector(0,m);	bigH=dvector(0,m);
	C1=dvector(0,m);	C2=dvector(0,m);
	C3=dvector(0,m);	C4=dvector(0,m);
	//for(i=0;i<ndata;i++) bigA[ndata-1][i]=bigA[i][ndata-1]=0.;
	//for(i=0;i<ndata;i++) bigB[ndata-1][i]=bigB[i][ndata-1]=0.;
	//for(i=0;i<ndata;i++) big[ndata-1][i]=big[i][ndata-1]=0.;
	//for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) big[i][j]=0.;
	for(l=0;l<m;l++) { i=idi[l];	j=idj[l];
		big[l]=bigA[l]=bigB[l]=0.;
		bigB[l]=(z[j+1]-z[j])*mup[l]-(z[i+1]-z[i])*mu[l] +z[j]-z[i];
		bigA[l]=(z[j+1]-z[j])*(nup[l]-mup[l])
				-(z[i+1]-z[i])*(nu[l]-mu[l]);
		for(k=0;k<=norder;k++)
			big[l] +=bc[k]*( gsl_pow_int(bigA[l],k)
				*gsl_pow_int(bigB[l],norder-k) )/(k+1);
		big[l] *= mask[l];
	}

	//  Special treatment for the i=j contributions.
	/*
	for(i=0;i<ndata;i++) big[i][i]=mask[i][i]*gsl_pow_int(a[i]*tau,norder);
	for(i=0;i<ndata;i++) big[i][i]=big[i][i]*mask0[i][i];
	*/
	
	//  Compute numerator and denominator totals.
	sf_bot=sf_top=0.;
	for(k=0;k<m;k++) {
		sf_bot += mask[k];
		sf_top += big[k];
	}
	
	//  Compute uncertainties using partial derivaties wrt observables.
	for(l=0;l<m;l++) { i=idi[l];	j=idj[l];
		bigG[l]=bigH[l]=C1[l]=C2[l]=C3[l]=C4[l]=0.;
		for(k=0;k<norder;k++) bigH[l]+=
			((gsl_pow_int(bigA[l],k)*gsl_pow_int(bigB[l],norder-k-1)
				)*bc[k]*(norder-k))/(k+1);
		for(k=1;k<=norder;k++) bigG[l]+=
			((gsl_pow_int(bigA[l],k-1)*gsl_pow_int(bigB[l],norder-k)
				)*bc[k]*k)/(k+1);
		C1[l]=((nup[l]-mup[l])*bigG[l] +mup[l]*bigH[l])*mask[l];
		C2[l]=((mup[l]-nup[l])*bigG[l] +(1.-mup[l])*bigH[l])*mask[l];
		C3[l]=((mu[l]-nu[l])*bigG[l] -mu[l]*bigH[l])*mask[l];
		C4[l]=((nu[l]-mu[l])*bigG[l] +(1.-mu[l])*bigH[l])*mask[l];
	}

	free_dvector(bigA,0,m); free_dvector(bigB,0,m);
	free_dvector(bigG,0,m); free_dvector(bigH,0,m);
	free_dvector(mu,0,m);	free_dvector(mup,0,m);
	free_dvector(nu,0,m);	free_dvector(nup,0,m);

	//  Compute partial derivatives wrt observables.
	dSdz=ldvector(0,ndata+1);
	for(k=0;k<=ndata;k++) dSdz[k]=0.;
	for(l=0;l<m;l++) {	i=idi[l];	j=idj[l];
		dSdz[j+1] += C1[l];	//C1[i][j];
		dSdz[j  ] += C2[l];	//C2[i][j];
		dSdz[i+1] += C3[l];	//C3[i][j];
		dSdz[i  ] += C4[l];	//C4[i][j];
	}
	for(k=0;k<=ndata;k++) if(!finite(dSdz[k])) dSdz[k]=0.;
	/*
	for(k=0;k<=ndata;k++) { dSdz[k]=0.;
		if(k<ndata) for(i=0;i<ndata;i++)
			dSdz[k] += C2[i][k] + C4[k][i];
		if(k>0) for(i=0;i<ndata;i++)
			dSdz[k] += C1[i][k-1] + C3[k-1][i];
		//if(!finite(dSdz[k])) dSdz[k]=0.;
	}
	*/
	for(sf_unc=0.,k=0;k<=ndata;k++) sf_unc+=(dSdz[k]*dSdz[k])*dz[k]*dz[k];
	sf_unc=sqrt(sf_unc);

	//printf("%e\t%Le\t%Le\n",tau,sf_top/sf_bot,sf_unc/sf_bot);
	printf("%e\t%Le\t%Le\t%Le\n",tau,sf_bot,sf_top,sf_unc);
	*ans1=sf_bot;	*ans2=sf_top;	*ans3=sf_unc;

	free_ldvector(dSdz,0,ndata+1);
	free_ivector(bc,0,norder+1);
	free_dvector(C1,0,m);	free_dvector(C2,0,m);
	free_dvector(C3,0,m);	free_dvector(C4,0,m);
	free_dvector(big,0,m);	free_dvector(mask,0,m);
	free_dvector(a,0,ndata);	free_dvector(b,0,ndata);
	free_dvector(t,0,ndata);	free_dvector(dt,0,ndata);
	free_dvector(z,0,ndata);	free_dvector(dz,0,ndata);
	free_ivector(idi,0,m);	free_ivector(idj,0,m);
}


/*****************************************************************************
 *   Another method for calculating the structure function.  Unused now.
 ****************************************************************************/
void	zzigzag(double tau, int norder, int ndata, 
	long double w_t[], long double w_dt[],
	long double w_z[], long double w_dz[],
	double	*ans1, double *ans2, double *ans3)
{
	long double	*dSdz,sf_unc,sf_top,sf_bot;;
	double	*a,*b,*t,*dt,*z,*dz,**tlo,**thi,**mask,**mu,**nu,**mup,**nup;
	double	**big,**bigA,**bigB,**bigG,**bigH,**C1,**C2,**C3,**C4;
	int	i,j,k,n1,**mask0;
	int	*bc;

	//  Abbreviate the structure function order; binomial coefficients.
	n1=norder+1;
	bc=ivector(0,norder+1);
	binco(norder,bc);
	//for(i=0;i<=norder;i++) printf("%i\t",bc[i]);	printf("\n");

	//  Copy the data into working arrays.
	t=dvector(0,ndata);		dt=dvector(0,ndata);
	z=dvector(0,ndata);		dz=dvector(0,ndata);
	for(i=0;i<ndata;i++) {
		t[i]=w_t[i];		dt[i]=w_dt[i];
		//z[i]=w_z[i];		dz[i]=w_dz[i];
		t[i]=dat_t[i];		dt[i]=dat_dt[i];
		z[i]=dat_z[i];		dz[i]=dat_dz[i];
	}

	//  Optional logarithmic modes.
	if(flag_logx) {
		for(i=0;i<ndata;i++) dt[i]=dt[i]/t[i];
		for(i=0;i<ndata;i++) t[i]=log(t[i]);
	}
	if(flag_logy) {
		for(i=0;i<ndata;i++) dz[i]=dz[i]/z[i];
		for(i=0;i<ndata;i++) z[i]=log(z[i]);
	}

	//  Compute the gradients of the line segments.
	a=dvector(0,ndata);		b=dvector(0,ndata);
	for(i=0;i<ndata-1;i++) {
		a[i]=(z[i+1]-z[i])/(t[i+1]+t[i]);
		b[i]=z[i]-a[i]*t[i];
	}	
	a[ndata-1]=a[ndata-2];		b[ndata-1]=b[ndata-2];
	//a[ndata-1]=b[ndata-1]=0.;

	//  Compute the overlap time intervals.
	tlo=dmatrix(0,ndata,0,ndata);
	thi=dmatrix(0,ndata,0,ndata);
	mask=dmatrix(0,ndata,0,ndata);
	mask0=imatrix(0,ndata,0,ndata);
	for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) {
		tlo[i][j]=gsl_max(t[j]+tau,t[i]);
		thi[i][j]=gsl_min(t[j+1]+tau,t[i+1]);
		mask[i][j]=thi[i][j]-tlo[i][j];
		mask0[i][j]=(mask[i][j]>0.)?1:0;
		mask[i][j]=(mask[i][j]>0.)?mask[i][j]:0.;
	}

	//  Compute fractional intervals of the overlapping lines.
	mu=dmatrix(0,ndata,0,ndata); mup=dmatrix(0,ndata,0,ndata);
	nu=dmatrix(0,ndata,0,ndata); nup=dmatrix(0,ndata,0,ndata);
	for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) {
		mu[i][j]=nu[i][j]=mup[i][j]=nup[i][j]=0.;
		if(mask0[i][j]==0) continue;
		if(i<ndata-1) {
			mu[i][j]=(tlo[i][j]-t[i])/(t[i+1]-t[i]);
			nu[i][j]=(thi[i][j]-t[i])/(t[i+1]-t[i]);
		}
		if(j<ndata-1) {
			mup[i][j]=(tlo[i][j]-t[j]-tau)/(t[j+1]-t[j]);
			nup[i][j]=(thi[i][j]-t[j]-tau)/(t[j+1]-t[j]);
		}
	}
	/*
	printf("%e < mu  <%e\n",findmin(mu,i,i,mask0),findmax(mu,i,i,mask0));
	printf("%e < nu  <%e\n",findmin(nu,i,i,mask0),findmax(nu,i,i,mask0));
	printf("%e < mup <%e\n",findmin(mup,i,i,mask0),findmax(mup,i,i,mask0));
	printf("%e < nup <%e\n",findmin(nup,i,i,mask0),findmax(nup,i,i,mask0));
	*/
	free_dmatrix(tlo,0,ndata,0,ndata); free_dmatrix(thi,0,ndata,0,ndata);

	/*  Terms in the SF contributions from each line segment.  */
	big=dmatrix(0,ndata,0,ndata);
	bigA=dmatrix(0,ndata,0,ndata);	bigB=dmatrix(0,ndata,0,ndata);
	bigG=dmatrix(0,ndata,0,ndata);	bigH=dmatrix(0,ndata,0,ndata);
	C1=dmatrix(0,ndata,0,ndata);	C2=dmatrix(0,ndata,0,ndata);
	C3=dmatrix(0,ndata,0,ndata);	C4=dmatrix(0,ndata,0,ndata);
	for(i=0;i<ndata;i++) bigA[ndata-1][i]=bigA[i][ndata-1]=0.;
	for(i=0;i<ndata;i++) bigB[ndata-1][i]=bigB[i][ndata-1]=0.;
	for(i=0;i<ndata;i++) big[ndata-1][i]=big[i][ndata-1]=0.;
	for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) big[i][j]=0.;
	for(i=0;i<ndata-1;i++) for(j=0;j<ndata-1;j++) {
		big[i][j]=bigA[i][j]=bigB[i][j]=0.;
		if(mask0[i][j]==0) continue;
		bigB[i][j]=(z[j+1]-z[j])*mup[i][j]-(z[i+1]-z[i])*mu[i][j]
				+z[j]-z[i];
		bigA[i][j]=(z[j+1]-z[j])*(nup[i][j]-mup[i][j])
				-(z[i+1]-z[i])*(nu[i][j]-mu[i][j]);
		big[i][j]=0.;
		for(k=0;k<=norder;k++)
			big[i][j] +=bc[k]*( gsl_pow_int(bigA[i][j],k)
				*gsl_pow_int(bigB[i][j],norder-k) )/(k+1);
		big[i][j] *= mask[i][j];
	}

	//  Special treatment for the i=j contributions.
	for(i=0;i<ndata;i++)
	//	big[i][i]=mask[i][i]*gsl_pow_int(a[i]*tau,norder);
		big[i][i]=big[i][i]*mask0[i][i];
	
	//  Compute numerator and denominator totals.
	sf_bot=sf_top=0.;
	for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) if(mask0[i][j]!=0) {
		sf_bot += mask[i][j];
		sf_top += big[i][j];
	}
	
	//  Compute uncertainties using partial derivaties wrt observables.
	for(i=0;i<ndata;i++) for(j=0;j<ndata;j++) {
		bigG[i][j]=bigH[i][j]=C1[i][j]=C2[i][j]=C3[i][j]=C4[i][j]=0.;
		if(mask0[i][j]==0) continue;
		for(k=0;k<norder;k++) bigH[i][j]+=
			((gsl_pow_int(bigA[i][j],k)
				*gsl_pow_int(bigB[i][j],norder-k-1)
					)*bc[k]*(norder-k))/(k+1);
		for(k=1;k<=norder;k++) bigG[i][j]+=
			((gsl_pow_int(bigA[i][j],k-1)
				*gsl_pow_int(bigB[i][j],norder-k)
					)*bc[k]*k)/(k+1);
		C1[i][j]=((nup[i][j]-mup[i][j])*bigG[i][j]
				+mup[i][j]*bigH[i][j])*mask[i][j];
		C2[i][j]=((mup[i][j]-nup[i][j])*bigG[i][j]
				+(1.-mup[i][j])*bigH[i][j])*mask[i][j];
		C3[i][j]=((mu[i][j]-nu[i][j])*bigG[i][j]
				-mu[i][j]*bigH[i][j])*mask[i][j];
		C4[i][j]=((nu[i][j]-mu[i][j])*bigG[i][j]
				+(1.-mu[i][j])*bigH[i][j])*mask[i][j];
	}

	free_dmatrix(bigA,0,ndata,0,ndata); free_dmatrix(bigB,0,ndata,0,ndata);
	free_dmatrix(bigG,0,ndata,0,ndata); free_dmatrix(bigH,0,ndata,0,ndata);
	free_dmatrix(mu,0,ndata,0,ndata); free_dmatrix(mup,0,ndata,0,ndata);
	free_dmatrix(nu,0,ndata,0,ndata); free_dmatrix(nup,0,ndata,0,ndata);

	/*  Compute partial derivatives wrt observables.  */
	dSdz=ldvector(0,ndata+1);
	for(k=0;k<=ndata;k++) { dSdz[k]=0.;
		if(k<ndata) for(i=0;i<ndata;i++)
			dSdz[k] += C2[i][k] + C4[k][i];
		if(k>0) for(i=0;i<ndata;i++)
			dSdz[k] += C1[i][k-1] + C3[k-1][i];
		//if(!finite(dSdz[k])) dSdz[k]=0.;
	}
	for(sf_unc=0.,k=0;k<=ndata;k++) sf_unc+=(dSdz[k]*dSdz[k])*dz[k]*dz[k];
	//for(sf_unc=0.,k=0;k<=ndata;k++) sf_unc+=gsl_pow_2(dSdz[k]*dz[k]);
	sf_unc=sqrt(sf_unc);

	//printf("%e\t%Le\t%Le\n",tau,sf_top/sf_bot,sf_unc/sf_bot);
	printf("%e\t%Le\t%Le\t%Le\n",tau,sf_bot,sf_top,sf_unc);
	*ans1=sf_bot;	*ans2=sf_top;	*ans3=sf_unc;

	free_ldvector(dSdz,0,ndata+1);
	free_dmatrix(C1,0,ndata,0,ndata);
	free_dmatrix(C2,0,ndata,0,ndata);
	free_dmatrix(C3,0,ndata,0,ndata);
	free_dmatrix(C4,0,ndata,0,ndata);
	free_dmatrix(big,0,ndata,0,ndata);
	free_dmatrix(mask,0,ndata,0,ndata);
	free_imatrix(mask0,0,ndata,0,ndata);
	free_dvector(a,0,ndata);	free_dvector(b,0,ndata);
	free_dvector(t,0,ndata);	free_dvector(dt,0,ndata);
	free_dvector(z,0,ndata);	free_dvector(dz,0,ndata);
	free_ivector(bc,0,norder+1);
}


/*****************************************************************************
 *  Attempt to open a file, and wait until successful.  This is sometimes
 *  useful when running the code on intermittently unstable disks.  :-/
 ****************************************************************************/
FILE	*wfopen(char fn[], char way[])
{	FILE	*f;
	while ((f=fopen(fn,way))==NULL) { sleep(30); }	return(f);
}


/*****************************************************************************
 *   Save the fitted parameters in binary format.
 ****************************************************************************/
void	save_parameters(char *fname, double *p, int npar)
{	FILE	*f;	int	i,n;		double	z;
	f=wfopen(fname,"w");
	//n=fwrite(&npar,sizeof(int),1,f);
	for(i=0;i<npar;i++) { z=p[i];	n=fwrite(&z,sizeof(double),1,f); }
	fclose(f);
}


/*****************************************************************************
 *  Read data for model fitting.  Format: T, dT, Z, dZ.  Count the columns.
 ****************************************************************************/
void	load_data(char *fname)
{
	FILE	*f;
	char	buf[200];
	double	t,dt,z,dz;
	int	i=0,n=0,tunset1=0,tunset2=0;
	if((f=fopen(fname,"r"))==NULL) {
		printf("FILE %s UNREADABLE\n",fname); 	exit(0);
	}
	while(fgets(buf,200,f) != NULL) {
		if(buf[0]=='#') continue;
		sscanf(buf,"%lf %lf %lf %lf",&t,&dt,&z,&dz);
		if((t<tstart)||(t>tend)) continue;
		n++;
	}
	fclose(f);
	dat_t=ldvector(0,n);	dat_dt=ldvector(0,n);
	dat_z=ldvector(0,n);	dat_dz=ldvector(0,n);
	f=wfopen(fname,"r");
	i=0;
	tunset1=(tmin>= FLT_MAX);
	tunset2=(tmax<=-FLT_MAX);
	while(fgets(buf, 200, f) != NULL) {
		if(buf[0]=='#') continue;
		dz=0.;		t=DBL_MAX;
		sscanf(buf,"%lf %lf %lf %lf",&t,&dt,&z,&dz);
		if(tunset1) tmin=min(tmin,t-dt);	
		if(tunset2) tmax=max(tmax,t+dt);
		if((t<tstart)||(t>tend)) continue;
		//printf("%e\t%e\t%e\t%e\n",t,dt,z,dz);
		dat_t[i]=t;	dat_dt[i]=dt;
		dat_z[i]=z;	dat_dz[i]=dz;
		if(dz<=0.) dat_dz[i]=sqrt(z*dt*2.)/(2.*dt);
		//printf("%Le\t%Le\t%Le\t%Le\n"
		//,dat_t[i],dat_dt[i],dat_z[i],dat_dz[i]);
		i++;
	}
	ndata=i;
	fclose(f);
	printf("ndata = %i\n",ndata);
}


/*****************************************************************************
 *  Sort an array of sample tau values from low to high.
 ****************************************************************************/
void	d_sorter(double *tau, int n)
{
	gsl_vector_long_double *v =gsl_vector_long_double_alloc(n);
	gsl_permutation * p = gsl_permutation_alloc(n);
	long double	*xt;	xt=ldvector(0L,n);
	int	j;
	for(j=0;j<n;j++) gsl_vector_long_double_set(v,j,tau[j]);
	gsl_sort_vector_long_double_index(p,v);
	for(j=0;j<n;j++) xt[j]=tau[p->data[j]];
	for(j=0;j<n;j++) tau[j]=xt[j];
	free_ldvector(xt,0L,n);
	gsl_permutation_free(p);	gsl_vector_long_double_free(v);
}


/*****************************************************************************
 *  Main program.
 ****************************************************************************/
int	main(int argc, char *argv[])
{
	FILE	*f;
	int	i;

	//  Initialise random number generator, seeded from the clock.
	gsl_rng_env_setup(); T_rng=gsl_rng_default; r_rng=gsl_rng_alloc(T_rng);
	gsl_rng_set(r_rng,tnow());

	//  Default global settings.
	for(i=0;i<argc;i++) { strcat(wholecmd,argv[i]); strcat(wholecmd," "); }
	gethostname(hostname,(size_t)128);
	printf("HOST = %s\n",hostname);
	sprintf(profname,"%s.dat",hostname);

	//  Parse command-line options.  They override defaults.
	read_flags(argc,argv);		setup();

	//  Assign arrays for the results.
	sf_bot=dvector(0,ntries);	sf_top=dvector(0,ntries);
	sf_unc=dvector(0,ntries);	tau=dvector(0,ntries);

	//  Load the observational data.
	printf("LOAD: %s\n",argv[1]);
	load_data(argv[1]);

	//if(strstr(argv[0],"sf_zigzag")==NULL) exit(0);

	//  Choose the sequence of timescale, tau.
	for(i=0;i<ndata;i++) {
		if(tmin>= FLT_MAX) tmin=min(tmin,dat_t[i]-dat_dt[i]);
		if(tmax<=-FLT_MAX) tmax=max(tmax,dat_t[i]+dat_dt[i]);
	}
	printf("%e < t < %e\n",tmin,tmax);
	for(i=0;i<ntries/2;i++) tau[i]=tmin*pow(tmax/tmin,i/(ntries/2-1.));
	for(i=ntries/2;i<ntries;i++)
		tau[i]=tmin+(tmax-tmin)*(i-ntries/2+1)/(ntries-ntries/2+1);
	d_sorter(tau,ntries);

	//  Execute the actual time series calculations.
	if(flag_pds) {
	//  Compute Leahy normalised power spectra at desired tau values.
		for(i=0;i<ntries;i++)
		zigpds(tau[i],sfn,ndata,dat_t,dat_dt,dat_z,dat_dz
				,&sf_bot[i],&sf_top[i],&sf_unc[i]);
	} else {
	//  Calculate structure functions at the desired tau values.
		for(i=0;i<ntries;i++)
		zigzag(tau[i],sfn,ndata,dat_t,dat_dt,dat_z,dat_dz
				,&sf_bot[i],&sf_top[i],&sf_unc[i]);
	}

	//  Save the results:
	//
	//  	column 1	=  tau, time-scale
	//  	column 2	=  W, denominator of SF
	//  	column 3	=  W*S, numerator of SF
	//  	column 4	=  W*dS, uncertainty in (W S)
	//
	//  	structure function (tau) = ( W*S +/- W*dS) / W
	//
	f=wfopen(profname,"w");
	fprintf(f,"#%s\n",wholecmd);
	for(i=0;i<ntries;i++) if(finite(sf_top[i]+sf_unc[i]))
		fprintf(f,"%e\t%e\t%e\t%e\n",tau[i]
			,sf_bot[i],sf_top[i],sf_unc[i]);
	fclose(f);
	
	//  Deallocate the data and output array variables.
	free_dvector(sf_bot,0,ntries);	free_dvector(sf_top,0,ntries);
	free_dvector(sf_unc,0,ntries);	free_dvector(tau,0,ntries);
	free_ldvector(dat_t,0,ndata);	free_ldvector(dat_dt,0,ndata);
	free_ldvector(dat_z,0,ndata);	free_ldvector(dat_dz,0,ndata);
	return(0);
}

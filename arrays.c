void	fatal(char msg[])
{	fprintf(stderr,"FATAL ERROR:\n%s\n",msg);	exit(1); }

int	*ivector(long n0, long N)
{	int	*V;
	V=(int *)malloc((size_t) ((N-n0+2)*sizeof(int)));
	if(!V) fatal("indigestion in ivector()"); else return(V-n0+1);
}

double	*dvector(long n0, long N)
{	double	*V;
	V=(double *)malloc((size_t) ((N-n0+2)*sizeof(double)));
	if(!V) fatal("indigestion in dvector()"); else return(V-n0+1);
}

long double	*ldvector(long n0, long N)
{	long double	*V;
	V=(long double *)malloc((size_t) ((N-n0+2)*sizeof(long double)));
	if(!V) fatal("indigestion in ldvector()"); else return(V-n0+1);
}

double	**dmatrix(long i0, long NI, long j0, long NJ)
{	double	**M;
	long	i,nx=NI-i0+1,ny=NJ-j0+1;
	M=(double **)malloc((size_t)((nx+1)*sizeof(double*)));
	if(!M) fatal("indigestion I in dmatrix()");
	M+=1;		M-=i0;
	M[i0]=(double *)malloc((size_t)((nx*ny+1)*sizeof(double)));
	if(!M[i0]) fatal("indigestion J in dmatrix()");
	M[i0]+=1;	M[i0]-=j0;
	for(i=i0+1;i<=NI;i++) M[i]=M[i-1]+ny;
	return(M);
}

int	**imatrix(long i0, long NI, long j0, long NJ)
{	int	**M;
	long	i,nx=NI-i0+1,ny=NJ-j0+1;
	M=(int **)malloc((size_t)((nx+1)*sizeof(int*)));
	if(!M) fatal("indigestion I in imatrix()");
	M+=1;		M-=i0;
	M[i0]=(int *)malloc((size_t)((nx*ny+1)*sizeof(int)));
	if(!M[i0]) fatal("indigestion J in imatrix()");
	M[i0]+=1;	M[i0]-=j0;
	for(i=i0+1;i<=NI;i++) M[i]=M[i-1]+ny;
	return(M);
}

void	free_ivector(int *V, long n0, long N)
{	free((char*)(V+n0-1)); }

void	free_dvector(double *V, long n0, long N)
{	free((char*)(V+n0-1)); }

void	free_ldvector(long double *V, long n0, long N)
{	free((char*)(V+n0-1)); }

void free_dmatrix(double **M, long i0, long NI, long j0, long NJ)
{	free((char*)(M[i0]+j0-1));	free((char*)(M+i0-1)); }

void free_ldmatrix(long double **M, long i0, long NI, long j0, long NJ)
{	free((char*)(M[i0]+j0-1));	free((char*)(M+i0-1)); }

void free_imatrix(int **M, long i0, long NI, long j0, long NJ)
{	free((char*)(M[i0]+j0-1));	free((char*)(M+i0-1)); }

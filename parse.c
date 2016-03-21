/*****************************************************************************
 *  Parse the command-line arguments as flag options, numerical parameters
 *  or equations to overwrite named variables.
 *****************************************************************************/
/*
	-mute	suppress some screen messages
	-tmp	put temporary files in /tmp instead of /unsafe
	-NOW	soften the imposition of curfew
	-rel	minimise max relative error rather than max absolute error
	-logx	take logarithms of X data, and adjust uncertainties
	-logy	take logarithms of Y data, and adjust uncertainties
	-pds	generate power spectrum with Leahy normalisation, not a SF.

	n	= order of structure function
	start	= omit data before this starting time
	end	= omit data after this ending time
	tmin	= minimum timescale for SF or PDS calculation
	tmax	= maximum timescale for SF or PDS calculation
	nx	= number of coefficients in X independent variable
	ny	= number of coefficients in Y independent variable
	guess	= maximum initial guess for free parameters
	try	= number of fitting attempts, before saving the best
*/

double	numin[42];	/*  Numerical command-line options.  */
long double	eqnrhs[42];	/*  Equation command-line options, values.  */
char	eqnlhs[42][80];	/*  Equation command-line options, variable labels.  */
int	numpars;
int	numequations=0;

int	is_number(char *test)
{
/*  Test whether a string represents a valid number.  */
	int c,bad,point,e;
	bad=point=e=0;
	for(c=0;test[c]!='\0';c++) {
		if( (!isdigit(test[c])) && (test[c]!='.')
		   && (test[c]!='e') &&(test[c]!='E')
		   && (test[c]!='+') &&(test[c]!='-') ) bad++;
		if(test[c]=='.') point++;
		if((test[c]=='e')||(test[c]=='E')) e++;
	}
	if((bad>0)||(point>1)||(e>1)) return(FALSE);
	else return(TRUE);
}


/***************************************************************************/

void	read_flags(argc,argv)
int	argc;
char	*argv[];
{
	char	par[80],*val,*wo,*end,*endp;
	long double	nval;
	int	p,nn,nf,neqn,l;
	numpars=0;
	neqn=nn=nf=0;
	for(p=1;p<argc;p++) {
		val=strchr(argv[p],'=');
		if(val!=NULL) {
			/*  This option is an equation.  */
			printf("equation:\t");
			eqnrhs[neqn]=nval=atof(val+1);
			l=val-argv[p];
			strcpy(par,argv[p]);
			val=strchr(par,'=');
			*val='\0';
			strcpy(eqnlhs[neqn],par);
			printf("%s = %Lf\n",par,nval);
			neqn++;
			continue;
		} else {
		if (is_number(argv[p])) {
			/*  This option is a number.  */
			numin[nn]=atof(argv[p]);
			printf("number : %e\n",numin[nn]);
			nn++;
		} else {
			/*  This option is a string or flag.  */
			l=strlen(argv[p]);
			/*
			*/
			if(argv[p][0]=='^') {	// file name
				int i;
				for(i=1;i<l;i++) profname[i-1]=argv[p][i];
				profname[i-1]='\0';
				printf("file : \"%s\"\n",profname);
				continue;
			} else
				printf("string : \"%s\"\t%i long\n",argv[p],l);

			if(!strcmp(argv[p],"-rel"))	flag_rel=TRUE;
			if(!strcmp(argv[p],"-mute"))	mute_setup=TRUE;
			if(!strcmp(argv[p],"-pds"))	flag_pds=TRUE;
			if(!strcmp(argv[p],"-logx"))	flag_logx=TRUE;
			if(!strcmp(argv[p],"-logy"))	flag_logy=TRUE;

			nf++;
		}
		}
	}
	numpars=nn;
	numequations=neqn;
	if(argc>1) {
		printf("===== program options =====\n");
		printf("\tequations : %i\n",numequations);
		printf("\tnumerical : %i\n",numpars);
		printf("\tcharacter : %i\n",argc-numpars-1);
	}
	return;
}


/***************************************************************************/

void	read_indat()
{
/*  Append equations stated in the "indat" setup file.  */
	FILE *indat;
	int	i;
	char	s[40],olds[40];
	long double	x;
	
	sprintf(olds,"ZYXWV");
	indat=fopen("indat","r");

	do {
		strcpy(olds,s);
		i=fscanf(indat,"%s = %Lf",s,&x);
		if(s[0]==EOF) break;
		if(!strcmp(s,".")) break;
		if(!mute_setup) printf("\t\"%s\"\t= % Le\n",s,x);
		strcpy(eqnlhs[numequations],s);
		eqnrhs[numequations]=x;
		numequations++;
	} while(strcmp(s,olds));

	fclose(indat);
}


/****************************************************************************/

void	setup()
{
	void show(char *s, long double x) { printf("\t\"%s\"\t= % Lg\n",s,x); }
	long double	x;
	char	s[40];
	int	i,p=0;
	for(i=numequations-1;i>=0;i--) {
		sprintf(s,"%s",eqnlhs[i]);
		x=eqnrhs[i];
		if(!strcmp(s,"nx" )) { nx=x; show(s,x); }
		if(!strcmp(s,"ny" )) { ny=x; show(s,x); }
		if(!strcmp(s,"guess" )) { aguess=x; show(s,x); }
		if(!strcmp(s,"try" )) { ntries=x; show(s,x); }
		if(!strcmp(s,"start" )) { tstart=x; show(s,x); }
		if(!strcmp(s,"end" )) { tend=x; show(s,x); }
		if(!strcmp(s,"tmin" )) { tmin=x; show(s,x); }
		if(!strcmp(s,"tmax" )) { tmax=x; show(s,x); }
		if(!strcmp(s,"n" )) { sfn=x; show(s,x); }
	}
}

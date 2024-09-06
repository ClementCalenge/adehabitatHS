#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void multpoco(double **tab, double *poco);
void aleadistrivec(double *vec, double *no);
void randksel(double *pu, int *nani, int *ni);
void rks(int *fac, double *pdsu, int *nani, int *nbani, int *nl);
void ksel(double *tab, int *fac, double *poidsut, int *nhab, 
	  int *nani, int *nloctot, double *ut, double *di,
	  double *marg, int *nombreani, double *eigenvp, 
	  double *poidsco, int *ewa);
void permutksel(double *tab, int *fac, double *poidsut, int *nhab,
		int *nani, int *nloctot, double *ut, double *di,
		double *marg, int *nombreani, int *npermut, 
		double *obseig, double *simeig, double *obsmarg,
		double *simmarg, double *eigenvp, 
		double *simtout, double *poidsco,
		int *ewa);
void calcsim(double *pix, double **pts, double *rg, 
	     int *nvar, int *npts, double *similarite);
void fctdomain(double *kascmod, double *ptsmod, 
	       double *range, int *npts, int *npix,  
	       int *nvar, double *qualhab);
void wml(double **used, double **avail, double *wmla, int na, int nh,
	 double **proj1, double **proj2, double *nbassocie, int krep);
void aclambda(double *util, double *dispo, int *nani, int *nhab,  
	      double *xxtxxtmod1, double *xxtxxtmod2, double *rnv,
	      double *wmla, int *nrep, double *wm, double *nb);
void rankma(double *used, double *avail, double *rankmap, double *rankmam,
	    double *rankmav, double *rankmanb, 
	    int *nhab, int *nani, int *nrep, double *rnv);
void matmudemi(double **X, double **Y);
void matmudemir(double *Xr, double *Yr, int *ncr);
void enfa(double **Z, double *p, int *nvar, int *npix,
	  double *vp);
void enfar(double *Zr, double *pr, int *nvar, int *npix,
	   double *vpr);
void randenfa(double **Z, double *p, int *nrep, double *res);
void randenfar(double *Zr, double *pr, int *nvar, int *npix,
	       int *nrep, double *resr);
void engen2008Ir(double *avr, double *usr, int *nliga, int *nligu, 
		 int *ncol, double *resr, int *nsimrar);
void engen2008r(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, int *idr, int *nidr, int *nsimr, 
		double *resr, int *nsimrar);
void free_ivector(int *v, int nl);
int invers(double **a, int n, double **b, int m);






/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    GetRNGstate();
    w = unif_rand();
    PutRNGstate();
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    for (i=1;i<=n;i++) {
	GetRNGstate();
	alea[i] = ((int) (1e8)*(unif_rand()));
	PutRNGstate();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.00000000001) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutees par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i;
    
    *tab = R_Calloc(l1+1, double *);
    for (i=0;i<=l1;i++) {
	*(*tab+i)= R_Calloc(c1+1, double);
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    *vec = R_Calloc(n+1, double);
    **vec = n;
    return;
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    *vec = R_Calloc(n+1, int);
    **vec = n;
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	R_Free(*(tab+i));
    }
    R_Free(tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    R_Free(vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    R_Free(vec);
    
}













/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitat               *****
 *********               -------------------------               *****
 *********************************************************************
 *********************************************************************
 */




/* ------------------------------------------------------ */
/* Multiply a table by the square root of columns weights */

void multpoco(double **tab, double *poco)
{

    /* Declaration of local variables */
    int nc, nl, i, j;
    double k;
    
    nl = tab[0][0];
    nc = tab[1][0];
    
    /* Main operations */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    k = poco[j];
	    tab[i][j] = tab[i][j]*sqrt(k);
	}
    }
    
}


/* ------------------------------------------------------ */
/* Distributes "no" points randomly in a vector with p elements */

void aleadistrivec(double *vec, double *no)
{
    
/* Declaration of local variables */
    double tmp, i, j, n, lv;
    
    n = *no;
    lv = vec[0];
    
    /* The random distribution of points in the elements of the vector */
    for (i=1; i<=n; i++) {
	tmp = alea();
	for (j=1; j<=lv; j++) {
	    if ((tmp >= (j-1)/lv)&&(tmp < j/lv))
		vec[(int) j]++;
	}
    }
}





/* *****************************************************
   
Randomly distributes the points within the home range 
of animals (used for K-select analysis)

***************************************************** */

void randksel(double *pu, int *nani, int *ni)
{
    /* Declaration of local variables */
    int i, j, k, l;
    double su, *tmp;
    
    l=1;
    
    /* For each animal, randomizes the points */
    for (k=1; k<=*nani; k++) {
	
	vecalloc(&tmp, ni[k]);
	su = 0;
	
	for (i=1; i<=ni[k];i++) {
	    su = su+pu[l];
	    l++;
	}
	
	aleadistrivec(tmp, &su);
	
	j=1;
	for (i=(l-ni[k]); i<l; i++) {
	    pu[i]=tmp[j];
	    j++;
	}
	
	freevec(tmp);
    }
    
}



/* *****************************************************

        The version to test randksel with R

   ***************************************************** */

void rks(int *fac, double *pdsu, int *nani, int *nbani, int *nl)
{

    /* Declaration and copy of R objects in local variables */
    
    int i, *fa, *ni;
    double *pu;
    
    vecalloc(&pu, *nl);
    vecintalloc(&fa, *nl);
    vecintalloc(&ni, *nani);
    
    for (i=1; i<=*nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=*nl; i++) {
	pu[i] = pdsu[i-1];
    }
    
    for (i=1; i<=*nani; i++) {
	ni[i] = nbani[i-1];
    }
    
    /* The function randksel */
    randksel(pu, nani, ni);
    
    /* Output */
    for (i=1; i<=*nl; i++) {
	pdsu[i-1] = pu[i];
    }
    
    /* Free memory */
    freevec(pu);
    freeintvec(fa);
    freeintvec(ni);
}



/* ****************************************************************
   *                                                              *
   *              K-select Analysis                               *
   *                                                              *
   **************************************************************** */


void ksel(double *tab, int *fac, double *poidsut, int *nhab, 
	  int *nani, int *nloctot, double *ut, double *di,
	  double *marg, int *nombreani, double *eigenvp, 
	  double *poidsco, int *ewa)
{
    
    /* Declaration of local variables */
    int i,j,k, sommeloctot;
    double **ta, *pu, **use, **ava, **mar, *poidsani;
    double **inertie, *valpro, *spu, *poco;
    int nh, na, nl, rang;
    int *fa, *ni;
    
    
    /* Memory Allocation for local variables */
    nl = *nloctot; /* total number of pixels */
    na = *nani; /* number of animals */
    nh = *nhab; /* number of variables */
    
    vecintalloc (&fa, nl); /* factor with one level per animal */
    vecalloc (&pu, nl); /* utilization weight of the pixels */
    vecalloc (&poidsani, na); /* weight of each animal (proportional
				 to the nb of locs) */
    vecalloc (&valpro, nh); /* eigenvalues of the analysis */
    vecalloc (&spu, na); /* utilization weights (sum equal to 1) */
    vecalloc(&poco, nh); /* column weights in the analysis (habitat weights) */
    taballoc (&ta, nl, nh); /* initial table */
    taballoc (&use, na, nh); /* table of used means */
    taballoc (&ava, na, nh); /* table of available means */
    taballoc (&mar, na, nh); /* marginality table */
    taballoc (&inertie, nh, nh); /* Inertia matrix */
    vecintalloc(&ni, na); /* number of available pixels per animal */
    

    
    /* Copy of the R objects into the local C variables */
    k = 0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nh; j++) {
	    ta[i][j] = tab[k];
	    k = k + 1;
	}
    }
    
    for (i=1; i<=nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=nl; i++) {
	pu[i] = poidsut[i-1];
    }
    
    for (i=1; i<=nh; i++) {
	poco[i] = poidsco[i-1];
    }
    
    for (i=1; i<=na; i++) {
	spu[i] = 0;
    }
    
    for (i=1; i<=na; i++) {
	ni[i] = nombreani[i-1];
    }
    
    
    /* Computes the number of relocs per animal */
    for (i=1; i<=na; i++) {
	for (k=1; k<=nl; k++) {
	    if (fa[k]==i) {
		spu[i] = spu[i] + pu[k];	
	    }
	}
    }
    
    
    /* Computes the total number of relocations */
    sommeloctot=0;
    for (i=1; i<=na; i++) {
	sommeloctot = sommeloctot + spu[i];	
    }
    
    
    /* Computes the used mean and available mean */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    for (k=1; k<=nl; k++) {
		if (fa[k]==i) {
		    ava[i][j]= ava[i][j] + ta[k][j]/ni[i];
		    use[i][j]=use[i][j] + (ta[k][j] * pu[k]);
		}
	    }
	}
    }
    
    /* and the marginality */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    use[i][j] = use[i][j] / spu[i];
	    mar[i][j] = use[i][j] - ava[i][j];
	}
    }
    
    /* Output toward the R object */
    k = 0;
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    ut[k] = use[i][j];
	    di[k] = ava[i][j];
	    marg[k] = mar[i][j];
	    k = k + 1;
	}
    }
    
    /* column weighting */
    multpoco(mar, poco);
    
    
    /* Computation of the weight of each animal */
    for (i=1; i<=na; i++) {
	if (*ewa==0)
	    poidsani[i] = (double) spu[i] / sommeloctot;
	if(*ewa==1)
	    poidsani[i] = (double) 1 / na;
    }
    
    sqrvec(poidsani);
    
    
    /* Computation of the eigenstructure */
    for (i=1; i<=na; i++) {
	for (j=1; j<=nh; j++) {
	    mar[i][j] = mar[i][j] * poidsani[i];
	}
    }
    prodmatAtAB(mar, inertie);
    DiagobgComp(nh, inertie, valpro, &rang);
    
    
    /* Output toward R */
    for (i = 1; i<=rang; i++) {
	eigenvp[i-1] = valpro[i];
    }
    
    
    /* Free memory */
    freeintvec (fa);
    freeintvec (ni);
    freevec (pu);
    freevec (poidsani);
    freevec(valpro);
    freevec (spu);
    freevec(poco);
    freetab (ta);
    freetab (use);
    freetab (ava);
    freetab (mar);
    freetab (inertie);
}




/* ****************************************************************
   *                                                              *
   *         Randomisation test related to the K-select           *
   *                                                              *
   **************************************************************** */


void permutksel(double *tab, int *fac, double *poidsut, int *nhab,
		int *nani, int *nloctot, double *ut, double *di,
		double *marg, int *nombreani, int *npermut, 
		double *obseig, double *simeig, double *obsmarg,
		double *simmarg, double *eigenvp, double *simtout, 
		double *poidsco, int *ewa)
{
    /* Declaration of local variables */
    double **ta, *pu, *obstout;
    int na, nh, nl, i, j, k, l, m, q, *ni, *fa, *numero, nbperm;
    
    /* Memory Allocation */
    na=*nani;
    nh=*nhab;
    nl=*nloctot;
    nbperm=*npermut;
    
    taballoc(&ta, nl, nh);
    vecintalloc(&fa, nl);
    vecintalloc(&numero, nbperm);
    vecalloc(&pu, nl);
    vecalloc(&obstout, na*nh);
    vecintalloc(&ni, na);
    
    
    /* One copies R objects toward C objects (local variables) */
    k = 0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nh; j++) {
	    ta[i][j] = tab[k];
	    k = k + 1;
	}
    }
    
    for (i=1; i<=nl; i++) {
	fa[i] = fac[i-1];
    }
    
    for (i=1; i<=nl; i++) {
	pu[i] = poidsut[i-1];
    }
    
    for (i=1; i<=na; i++) {
	ni[i] = nombreani[i-1];
    }
    
    /* --------------------------------------- */
    /* Main Computations */
    
    /* Observed kselect */
    ksel(tab,  fac, poidsut, nhab, 
	 nani, nloctot, ut, di,
	 marg, nombreani, eigenvp, poidsco, ewa);
    
    
    /* observed values are placed in the output */
    k=0;
    j=0;
    *obseig = eigenvp[0];
    for (i=0; i<(na*nh); i++) {
	obsmarg[j] = obsmarg[j] + (marg[i] * marg[i] * poidsco[k]);
	if (k==(nh-1)) {
	    k=-1;
	    j++;
	}
	k++;
    }
    
  
    /* Coordinates of the marginality vectors */
    for (i=1; i<=na*nh; i++) {
	obstout[i] = marg[i-1];
    }
    
    
    
    /* The permutations */
    m=0; /* will indicate the step for simmarg */
    q=0; /* will indicate the step for simtout */
    
    for (k=1; k<=nbperm; k++) {
	
	/* One permutes */
	randksel(pu, nani, ni);
	
	/* One copies the randomized utilization weights in poidsut */
	for (i=0; i<nl; i++) {
	    poidsut[i]=pu[i+1];
	}
	
	/* and therefore, go */
	ksel(tab,  fac, poidsut, nhab, 
	     nani, nloctot, ut, di,
	     marg, nombreani, eigenvp, poidsco, ewa);
	
	
	/* Simulated values in the output */
	simeig[k-1] = eigenvp[0];
	l=0;
	
	for (i=0; i<(na*nh); i++) {
	    simmarg[m] = simmarg[m] + (marg[i] * marg[i] * poidsco[l]);
	    if (l==(nh-1)) {
		l=-1;
		m++;
	    }
	    l++;
	}
	
	for (i=0; i<(na*nh); i++) {
	    simtout[q] = marg[i];
	    q++;
	}
    }
    
    /* obstout is put again in marg */
    for (i=0; i<(na*nh); i++) {
	marg[i] = obstout[i+1];
    }
    
    
    /* Memory free */
    freetab(ta);
    freeintvec(fa);
    freevec(pu);
    freevec(obstout);
    freeintvec(ni);
    freeintvec(numero);
}









/* ****************************************************************
   *                                                              *
   *   DOMAIN: estimation of the potential distribution range     *
   *                                                              *
   **************************************************************** */

void calcsim(double *pix, double **pts, double *rg, 
	     int *nvar, int *npts, double *similarite)
{
    /* Declarations of variables */
    int no,nv, i, j;
    double *vecqual, *temp, nib;
    
    no = *npts;
    nv = *nvar;
    
    vecalloc(&vecqual, no);
    vecalloc(&temp, nv);
    
    /* Computation of the similarity */
    for (i=1; i<=no; i++) {
	nib = 0;
	for (j=1; j<=nv; j++) {
	    temp[j] = fabs(pix[j]-pts[i][j])/rg[j];
	    nib = nib + temp[j];
	}
	vecqual[i] = 1 - (1/((double) nv))*nib;
    }
    
    /* Computation of the habitat quality
       using the max of the similarity */
    *similarite = vecqual[1];
    
    for (i=2; i<=no; i++) {
	if (vecqual[i]>*similarite)
	    *similarite = vecqual[i];
    }
    
    
    /* Free memory */
    freevec(vecqual);
    freevec(temp);
    
}


/* For interaction with R */

void fctdomain(double *kascmod, double *ptsmod, 
	       double *range, int *npts, int *npix,  
	       int *nvar, double *qualhab)
{
    /* Declarations of variables */
    int no,np,nv, i, j, k;
    double **kasc, **pts, *rg, *pix, sim;
    
    /* Memory allocation */
    no = *npts;
    np = *npix;
    nv = *nvar;
    
    taballoc(&kasc, np, nv);
    taballoc(&pts, no, nv);
    vecalloc(&rg, nv);
    vecalloc(&pix, nv);
    
    /* R to C*/
    k = 0;
    for (i = 1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    kasc[i][j] = kascmod[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= no; i++) {
	for (j = 1; j <= nv; j++) {
	    pts[i][j] = ptsmod[k];
	    k++;
	}
    }
    
    for (i=1; i<=nv; i++) {
	rg[i] = range[i-1];
    }
    
    
    /* The core of the function */
    for (i=1; i<=np; i++) {
	for (j=1; j<=nv; j++) {
	    pix[j] = kasc[i][j];
	}
	calcsim(pix, pts, rg, &nv, &no, &sim);
	qualhab[i-1] = sim;
    }
    
    /* Free memory */
    freetab(kasc);
    freetab(pts);
    freevec(pix);
    freevec(rg);
    
}





/***********************************************************************
 *********                 Compositional analysis                 ******
 ***********************************************************************/


/* weighted mean lambda: compo analysis with R */

void wml(double **used, double **avail, double *wmla, int na, int nh,
	 double **proj1, double **proj2, double *nbassocie, int krep)
{
    /* Declaration of variables */
    double **dlr, *moydlr, *nadlr, **dlrtmp, **mod1, **mod2, **res1, **res2;
    double **SCEres1, **SCEres2, *vp1, *vp2, det1, det2, *vecalea, *aleamu;
    int i, j, k, idcol, *vecindice, rg1, rg2;
    int jb;
    
    /* Memory allocation */
    taballoc(&dlr, na, (nh*(nh-1)));
    taballoc(&mod1, na, (nh-1));
    taballoc(&mod2, na, (nh-1));
    taballoc(&SCEres1, (nh-1), (nh-1));
    taballoc(&SCEres2, (nh-1), (nh-1));
    taballoc(&dlrtmp, na, (nh-1));
    taballoc(&res1, na, (nh-1));
    taballoc(&res2, na, (nh-1));
    vecintalloc(&vecindice, nh-1);
    vecalloc(&nadlr, nh -1);
    vecalloc(&moydlr, nh-1);
    vecalloc(&vp1, nh-1);
    vecalloc(&vp2, nh-1);
    vecalloc(&aleamu, 2);
    vecalloc(&vecalea, na);
    
    aleamu[1] = 1;
    aleamu[2] = -1;
    
    jb = 0;
    
    /* Random permutation for each animal */
    for (i = 1; i <= na; i++) {
	aleapermutvec(aleamu);
	vecalea[i] = aleamu[1];
    }
    
    /* When krep == 1 : First repetition of the process: 
       Computation of the observed lambda */
    if (krep == 1) {
	for (i = 1; i<=na; i++) {
	    vecalea[i] = 1;
	}
    }
    
    
    /* empty nbassocie */
    for (i = 1; i <= nh; i++)
	nbassocie[i] = 0;
    
    
    /* loop to fill the DLR */
    for (k = 1; k <= nh; k++) {
	i = 1;
	
	/* build the vectors of indices */
	for (j = 1; j <= nh; j++) {
	    if (j != k) {
		vecindice[i] = j;
		i++;
	    }
	}
	
	/* Set the mean and the number of non-missing values to 0 */
	for (j = 1; j <= (nh-1); j++) {
	    moydlr[j] = 0;
	    nadlr[j] = 0;
	}
	
	/* First fill of the DLR */
	for (j = 1; j <= (nh-1); j++) {
	    jb = vecindice[j];
	    idcol = (nh - 1) * (k - 1) + j;
	    for (i = 1; i <= na; i++) {
		if ((fabs(avail[i][jb]) > 0.000000001)&&(fabs(avail[i][k]) > 0.000000001)) {
		    dlr[i][idcol] = (log(used[i][jb] / used[i][k]) - 
				     log(avail[i][jb] / avail[i][k])) * vecalea[i];
		    
		    /* computes the mean */
		    moydlr[j] = moydlr[j] + dlr[i][idcol]; 
		    nadlr[j]++;
		}
	    }
	}
	
	for (j = 1; j <= (nh-1); j++) {
	    moydlr[j] = moydlr[j] / nadlr[j];
	}
	
	/* Second loop: replace missing values */
	for (j = 1; j <= (nh-1); j++) {
	    idcol = (nh - 1) * (k - 1) + j;
	    jb = vecindice[j];
	    for (i = 1; i <= na; i++) {
		if ( (fabs(avail[i][jb]) < 0.000000001)||(fabs(avail[i][k])< 0.000000001))
		    dlr[i][idcol] = moydlr[j];
	    }
	}
	
	/* extraction of DLRtmp */
	for (i = 1; i <= na; i++) {
	    for (j = 1; j <= (nh-1); j++) {
		idcol = (nh - 1) * (k - 1) + j;
		dlrtmp[i][j] = dlr[i][idcol];
	    }
	}
	
	/* Computes the models */
	prodmatABC(proj1, dlrtmp, mod1);
	prodmatABC(proj2, dlrtmp, mod2);
	
	/* Computes the residuals */
	for (i = 1; i <= na; i++) {
	    for (j = 1; j <= nh-1; j++) {
		res1[i][j] = dlrtmp[i][j] - mod1[i][j];
		res2[i][j] = dlrtmp[i][j] - mod2[i][j];
	    }
	}
	
	/* The sum of squares */
	prodmatAtAB(res1, SCEres1);
	prodmatAtAB(res2, SCEres2);
	
	/* The determinant */
	DiagobgComp(nh-1, SCEres1, vp1, &rg1);
	DiagobgComp(nh-1, SCEres2, vp2, &rg2);
	det1 = 1;
	det2 = 1;
	
	if (rg1 != (nh-1)) {
	    det1 = 0;
	} else {
	    for (i = 1; i <= rg1; i++) {
		det1 = det1 * vp1[i];
	    }
	}
	
	if (rg2 != nh-1) {
	    det2 = 0;
	} else {
	    for (i = 1; i <= rg2; i++) {
		det2 = det2 * vp2[i];
	    }
	}
	
	wmla[k] = det1 / det2;
	for (i = 1; i <= (nh-1); i++)
	    nbassocie[k] = nbassocie[k] + nadlr[i];
    }
    
    /* free memory */
    freetab(dlr);
    freetab(mod1);
    freetab(mod2);
    freetab(SCEres1);
    freetab(SCEres2);
    freetab(dlrtmp);
    freetab(res1);
    freetab(res2);
    freeintvec(vecindice);
    freevec(nadlr);
    freevec(moydlr);
    freevec(vp1);
    freevec(vp2);
    freevec(aleamu);
    freevec(vecalea);
    
}




/* aclambda allows the computation of lambda in compositional analysis */

void aclambda(double *util, double *dispo, int *nani, int *nhab,  
	      double *xxtxxtmod1, double *xxtxxtmod2, double *rnv,
	      double *wmla, int *nrep, double *wm, double *nb)
{
    /* Declarations of variables */
    int na, nh, i, j, k, nr;
    double **ut, **di, **proj1, **proj2, *lilamb, *linb, sumnb;
    
    /* Memory allocation */
    na = *nani;
    nr = *nrep;
    nh = *nhab;
    
    taballoc(&ut, na, nh);
    taballoc(&di, na, nh);
    taballoc(&proj1, na, na);
    taballoc(&proj2, na, na);
    vecalloc(&lilamb, nh);
    vecalloc(&linb, nh);
    
    
    /* From R to C */
    /* use */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    ut[i][j] = util[k];
	    if (fabs(ut[i][j]) < 0.000000001)
		ut[i][j] = *rnv;
	    k++;
	}
    }
    
    /* availability */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    di[i][j] = dispo[k];
	    k++;
	}
    }
    
    /* projector 1 */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= na; j++) {
	    proj1[i][j] = xxtxxtmod1[k];
	    k++;
	}
    }
    
    /* projector 2 */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= na; j++) {
	    proj2[i][j] = xxtxxtmod2[k];
	    k++;
	}
    }
    
    
    /* Beginning of the loop */
    for (k = 1; k <= nr; k++) {
	/* weighted mean lambda */
	wml(ut, di, lilamb, na, nh,
	    proj1, proj2, linb, k);
	sumnb = 0.0;
	for (i = 1; i <= nh; i++)
	    sumnb = sumnb + linb[i];
	for (i = 1; i <= nh; i++)
	    wmla[k-1] = wmla[k-1] + ((lilamb[i] * linb[i]) / sumnb);
	
	/* C to R */
	if (k == 1) {
	    for (i = 1; i <= nh; i++) {
		wm[i-1] = lilamb[i];
		nb[i-1] = linb[i];
	    }
	}
    }
    
    
    /* Free memory */
    freetab(ut);
    freetab(di);
    freetab(proj1);
    freetab(proj2);
    freevec(lilamb);
    freevec(linb);
    
}


/* The ranking matrix for compositional analysis */

void rankma(double *used, double *avail, double *rankmap, double *rankmam,
	    double *rankmav, double *rankmanb, int *nhab, 
	    int *nani, int *nrep, double *rnv)
{
    /* Declarations of variables */
    int i, j, k, nh, na, nr, r;
    double **u, **a, **rmp, **rmm, **rmv, **rmnb;
    double *dlrtmp, *vecalea, val, moy;
    double *aleamu, **tabani;
    
    /* Memory Allocation */
    nh = *nhab;
    na = *nani;
    nr = *nrep;
    r = 0;
    
    taballoc(&u, na, nh);
    taballoc(&a, na, nh);
    taballoc(&rmv, nh, nh);
    taballoc(&rmp, nh, nh);
    taballoc(&rmm, nh, nh);
    taballoc(&rmnb, nh, nh);
    vecalloc(&dlrtmp, na);
    vecalloc(&vecalea, nr);
    vecalloc(&aleamu, 2);
    taballoc(&tabani, nr, na);
    aleamu[1] = -1;
    aleamu[2] = 1;
    
    /* Fill the table */
    k = 0;
    for (i = 1; i <= na; i++) {
	for (j = 1; j <= nh; j++) {
	    u[i][j] = used[k];
	    a[i][j] = avail[k];
	    if (fabs(u[i][j]) < 0.000000001)
		u[i][j] = *rnv;
	    k++;
	}
    }
    
    /* Fill the table tabani */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= na; j++) {
	    aleapermutvec(aleamu);
	    tabani[i][j] = aleamu[1];
	}
    }
    for (i = 1; i<=na; i++) {
	tabani[1][i] = 1;
    }
    
    /* Beginning of the loop */
    for (k = 1; k <= nh; k++) {
	for (j = 1; j <= nh; j++) {
	    for (r = 1; r <= nr; r++) {
		moy = 0;
		/* Fills first the DLR per animal */
		for (i = 1; i <= na; i++) {
		    if ((fabs(a[i][j])> 0.000000001)&&(fabs(a[i][k]) > 0.000000001)) {
			dlrtmp[i] = (log(u[i][j]/u[i][k]) - log(a[i][j]/a[i][k])) * tabani[r][i];
			moy = moy + dlrtmp[i];
			if (r == 1)
			    rmnb[j][k]++;
		    }
		}
		
		/* Computes the mean */
		moy = moy / rmnb[j][k];
		if (r==1)
		    rmv[j][k] = moy;
		vecalea[r] = moy;
	    }
	    
	    /* Computes P */
	    val = rmv[j][k];
	    for (r = 1; r <= nr; r++) {
		if (val < vecalea[r])
		    rmm[j][k]++;
		if (val > vecalea[r])
		    rmp[j][k]++;
	    }
	}
    }
    
    /* C to R */
    k = 0;
    for (i=1; i<=nh; i++) {
	for (j=1; j<=nh; j++) {
	    rankmap[k] = rmp[i][j];
	    rankmam[k] = rmm[i][j];
	    rankmav[k] = rmv[i][j];
	    rankmanb[k] = rmnb[i][j];
	    k++;
	}
    }
    
    
    /* Free memory */  
    freetab(rmv);
    freetab(rmp);
    freetab(rmm);
    freetab(rmnb);
    freevec(dlrtmp);
    freetab(u);
    freetab(a);
    freevec(vecalea);
    freevec(aleamu);
    freetab(tabani);
}







/* *******************************************

   Transformation of a square symetric matrix 
   into a matrix at the power -1/2
   
   ******************************************* */

void matmudemi(double **X, double **Y)
{
    /* Declaration of the variables */
    int i, j, nc, rg;
    double **U, **L, *lambda, **Ubis, **Uter;
    
    /* Memory Allocation */
    nc = X[0][0];
    taballoc(&U, nc, nc);
    taballoc(&Ubis, nc, nc);
    taballoc(&Uter, nc, nc);
    taballoc(&L, nc, nc);
    vecalloc(&lambda, nc);
    
    /* Fill Xtmp */
    for (i = 1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    U[i][j] = X[i][j];
	}
    }
    
    /* Eigenstructure of X */
    DiagobgComp(nc, U, lambda, &rg);
    
    /* Matrix lambda -1/2 */
    for (i = 1; i<=nc; i++) {
	L[i][i] = 1 / sqrt(lambda[i]);
    }
    
    /* Result */
    prodmatABC(U, L, Ubis);  
    for (i = 1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    Uter[i][j] = U[j][i];
	}
    }
    prodmatABC(Ubis, Uter, Y);
    
    /* Free memory */
    freetab(U);
    freetab(Ubis);
    freetab(Uter);
    freetab(L);
    freevec(lambda);
}



/* The same, for R */

void matmudemir(double *Xr, double *Yr, int *ncr)
{
    /* Declaration of variables */
    int i, j, k, nc;
    double **X, **Y;

    /* Memory allocation */
    nc = *ncr;
    taballoc(&X, nc, nc);
    taballoc(&Y, nc, nc);
    
    /* R to C */
    k = 0;
    for (i=1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    X[i][j] = Xr[k];
	    k++;
	}
    }
    
    /* matmudemi */
    matmudemi(X, Y);
    
    /* C to R */
    k = 0;
    for (i=1; i <= nc; i++) {
	for (j = 1; j <= nc; j++) {
	    Yr[k] = Y[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(X);
    freetab(Y);
}





/* *******************************************

   Enfa

   ****************************************** */

void enfa(double **Z, double *p, int *nvar, int *npix,
	  double *vp)
{
    /* Declaration of local variables */
    double *m, *z, *y, **W, **Rs, **Rg, **Zbis, **Rsm12, norz;
    double **Wtmp, **H, **Iv, **yyt, **Ivmyyt, **Htmp;
    int i, j, nv, np, rg;
    
    /* Memory allocation */
    np = *npix;
    nv = *nvar;
    norz = 0;
    rg = 0;
    
    taballoc(&Zbis, np, nv);
    vecalloc(&m, nv);
    vecalloc(&z, nv);
    vecalloc(&y, nv);
    taballoc(&W, nv, nv);
    taballoc(&Iv, nv, nv);
    taballoc(&Ivmyyt, nv, nv);
    taballoc(&Htmp, nv, nv);
    taballoc(&yyt, nv, nv);
    taballoc(&H, nv, nv);
    taballoc(&Wtmp, nv, nv);
    taballoc(&Rg, nv, nv);
    taballoc(&Rs, nv, nv);
    taballoc(&Rsm12, nv, nv);
    
    /* Marginality */
    for (j = 1; j<=nv; j++) {
	for (i = 1; i <= np; i++){
	    m[j] = m[j] + p[i] * Z[i][j];
	}
    }
    
    /* Rs and Rg */
    for (i = 1; i<=np; i++) {
	for (j = 1; j <= nv; j++) {
	    Zbis[i][j] = Z[i][j] * sqrt(p[i]);
	}
    }
    prodmatAtAB(Zbis, Rs);
    for (i = 1; i<=np; i++) {
	for (j = 1; j <= nv; j++) {
	    Zbis[i][j] = Z[i][j] * sqrt((1/ ((double) np)));
	}
    }
    prodmatAtAB(Zbis, Rg);
    
    /* Rs^-1/2  */
    matmudemi(Rs,Rsm12);
    
    /* z */
    for (i = 1; i <= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    z[i] = z[i] + Rsm12[i][j] * m[j];
	}
    }
    
    /* norm of z */
    for (i = 1; i <= nv; i++) {
	norz = norz + (z[i] * z[i]);
    }
    norz = sqrt(norz);
    
    /* y */
    for (i = 1 ; i <= nv; i++) {
	y[i] = z[i] / norz;
    }
    
    /* W */
    prodmatABC(Rsm12, Rg, Wtmp);
    prodmatABC(Wtmp, Rsm12, W);
    
    
    /* **************************** */
    /* The large part: H            */
    /* **************************** */
    
    /* yyt */
    
    for (i = 1; i<= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    yyt[i][j] = y[i] * y[j];
	}
    }
    
    
    /* Iv */
    
    for (i = 1; i <= nv; i++) {
	Iv[i][i] = 1;
    }
    
    
    /* Ivmyyt */
    
    for (i = 1; i <= nv; i++) {
	for (j = 1; j <= nv; j++) {
	    Ivmyyt[i][j] = Iv[i][j] - yyt[i][j];
	}
    }
    
    
    /* And finally, H */
    
    prodmatABC(Ivmyyt, W, Htmp);
    prodmatABC(Htmp, Ivmyyt, H);
    
    
    /* Eigenstructure of H */
    DiagobgComp(nv, H, vp, &rg);
    
    
    /* Free memory */
    freevec(m);
    freevec(z);
    freevec(y);
    freetab(W);
    freetab(Iv);
    freetab(Ivmyyt);
    freetab(Htmp);
    freetab(yyt);
    freetab(H);
    freetab(Wtmp);
    freetab(Rg);
    freetab(Rs);
    freetab(Rsm12);
    freetab(Zbis);
}


/* *****************************************************

ENFA with R

***************************************************** */

void enfar(double *Zr, double *pr, int *nvar, int *npix,
	   double *vpr)
{
    /* Declaration of variables */
    int i, j, k, np, nv;
    double **Z, *p, *vp;
    
    /* Memory allocation */
    taballoc(&Z, *npix, *nvar);
    vecalloc(&p, *npix);
    vecalloc(&vp, *nvar);
    
    np = *npix;
    nv = *nvar;
    
    /* R to C */
    k = 0;
    for (i=1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    Z[i][j] = Zr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= np; i++) {
	p[i] = pr[i-1];
    }
    
    /* ENFA ...*/
    enfa(Z, p, &nv, &np, vp);
    
    
    /* ... C to R ... */
    for (i = 1; i <= nv; i++) {
	vpr[i-1] = vp[i];
    }
    
    /* ... Free memory */
    freetab(Z);
    freevec(p);
    freevec(vp);
    
}



/* ********************************************************
   
Randomization in the ENFA: test of the first eigenvalue of specialization

******************************************************** */

void randenfa(double **Z, double *p, int *nrep, double *res)
{
    /* Declaration of variables */
    int i, j, k, nv, np, ntot;
    double *psim, *vp;
    
    /* Memory Allocation */
    np = Z[0][0];
    nv = Z[1][0];
    ntot = 0;
    vecalloc(&psim, np);
    vecalloc(&vp, nv);
    
    /* Counts the total number of points */
    for (i = 1; i <= np; i++) {
	ntot = ntot + p[i];
    }
    
    /* Beginning of the randomization porocess */
    for (k = 1; k <= *nrep; k++) {
	
	/* empty vector psim */
	for (i = 1; i <= np; i++) {
	    psim[i] = 0;
	}
	
	/* randomization of locs in the vector psim */
	for (i = 1; i <= ntot; i++) {
	    j = (int) (np * alea());
	    psim[j]++;
	}
	
	/* vector of weight */ 
	for (i = 1; i <= np; i++) {
	    psim[i] = psim[i] / ((double) ntot);
	}
	
	/* ... and ENFA */
	enfa(Z, psim, &nv, &np, vp); 
	
	/* storage in res... */
	res[k] = vp[1];
	
	/* ... and end of the loop */
    }
    
    /* Free memory */
    freevec(psim);
    freevec(vp);
}


/* The same but for external call from R */

void randenfar(double *Zr, double *pr, int *nvar, int *npix,
	       int *nrep, double *resr)
{
    /* Declaration of local variables */
    int i, j, k, nv, np, nr;
    double **Z, *p, *res;
    
    /* Memory Allocation */
    np = *npix;
    nv = *nvar;
    nr = *nrep;
    taballoc(&Z, np, nv);
    vecalloc(&p, np);
    vecalloc(&res, nr);
    
    /* R to C */
    k = 0;
    for (i=1; i <= np; i++) {
	for (j = 1; j <= nv; j++) {
	    Z[i][j] = Zr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= np; i++) {
	p[i] = pr[i-1];
    }
    
    /* C Function */
    randenfa(Z, p, &nr, res); 
    
    /* C to R */
    for (i = 1; i <= nr; i++) {
	resr[i-1] = res[i];
    }
    
    /* free memory */
    freevec(p);
    freevec(res);
    freetab(Z);
    
}




/* *********************************************************************
 *                                                                     *
 *                 Method of Engen (2008)                              *
 *                                                                     *
 ***********************************************************************/


/* 
   sources of the package deal-1.2-30. 
   Code from Susanne Gammelgaard Bottcher <alma@math.aau.dk>, 
   Claus Dethlefsen <cld@rn.dk>. 
   Useful for matrix inverse!
*/


int *ivector(int nl, int nh)
{
   int *v;

   v=(int *) R_alloc((unsigned) (nh-nl+1)*sizeof(int),sizeof(int));
   if ( v == NULL ){
      error("memory allocation failure in ivector()"); return(NULL);
   }
   return v-nl;
}

void free_ivector(int *v, int nl) { free((char*) (v+nl)); }



int invers(double **a, int n, double **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i,icol=1,irow=1,j,k,l,ll;
   double big,dum,pivinv;

   if( (indxc = ivector(1,n)) == NULL){ return(-1); }
   if( (indxr = ivector(1,n)) == NULL){ return(-1); }
   if( (ipiv  = ivector(1,n)) == NULL){ return(-1); }
   for (j=1;j<=n;j++) ipiv[j]=0;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if (ipiv[j] != 1)
            for (k=1;k<=n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(a[j][k]) >= big) {
                     big=fabs(a[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1){
                  error("Invers: Singular Matrix-1");
                  return(-1);
               }
            }
      ++(ipiv[icol]);
      if (irow != icol) {
         for (l=1;l<=n;l++){
            double temp=a[irow][l]; a[irow][l]=a[icol][l]; a[icol][l]=temp;
         }
         for (l=1;l<=m;l++){
            double temp=b[irow][l]; b[irow][l]=b[icol][l]; b[icol][l]=temp;
         }
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0){
         error("Invers: Singular Matrix-2");
         return(-1);
      }
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
         if (ll != icol) {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
         }
   }
   for (l=n;l>=1;l--) {
      if (indxr[l] != indxc[l]){
         for (k=1;k<=n;k++){
            double temp     = a[k][indxr[l]];
            a[k][indxr[l]] = a[k][indxc[l]];
            a[k][indxc[l]] = temp;
         }
      }
   }
   return(0);
}


/* end of the sources of deal */

/* Main function */


void engen2008r(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, int *idr, int *nidr, int *nsimr, 
		double *resr, int *nsimrar)
{
    /* declaration of variables */
    double **av, **us, **nsco, *varR, tmp, res, **var, mu1, mu2, **Akmo;
    double **inv1, **zer, *a, sigkk, *Wk, *tmp2, m, s, **nscob, **nscoav;
    double *obs, **sammr, vartot, sig2, *Zi, *thetai, tmp3;
    double **mu;
    int nla, nlu, nc, *id, i, j, k, nid, nsim, e, r,l, *index;
    int *indexR, nsimra, b, *ni;
    
    /* memory allocation */
    nc=*ncol;
    nla=*nliga;
    nlu=*nligu;
    nsim=*nsimr;
    nid=*nidr;
    nsimra=*nsimrar;
    e=0;
    r=1;
    tmp=0.0;
    tmp3=0.0;
    b=1;

    /* used and available */
    taballoc(&av, nla, nc);
    taballoc(&us, nlu, nc);
    
    /* the id */
    vecintalloc(&id, nlu);
    
    /* the results */
    taballoc(&mu, nsimra, (nc * 5) );
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	vecalloc(&Wk, (nc-1));
	taballoc(&Akmo, nc-1, nc-1);
	vecalloc(&a, (nc-1));
	vecalloc(&tmp2, (nc-1));
	vecintalloc(&index, nc);
	taballoc(&var, nc, nc);
	taballoc(&inv1, nc-1, nc-1);
	taballoc(&zer, nc-1, nc-1);
    }
    taballoc(&nscob, nlu, nc);

    /* elements required for the calculations
       of the statistics */
    taballoc(&sammr, nsim, nid);
    vecintalloc(&ni, nid);
    vecalloc(&Zi, nid);
    vecalloc(&thetai, nid);

    /* elements required for used normal scores */
    taballoc(&nsco, nlu, nc);
    
    /* elements required for available normal scores */
    taballoc(&nscoav, nla, nc);

    /* elements used in the R API */
    varR = (double *) R_alloc(nla, sizeof(double));
    indexR = (int *) R_alloc(nla, sizeof(int));



        
    /* R -> C objects */
    k=0;
    for (i=1; i<=nla; i++) {
	for (j=1; j<=nc; j++) {
	    av[i][j]=avr[k];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    us[i][j]=usr[k];
	    k++;
	}
    }
    
    for (i=1; i<=nlu; i++) {
	id[i]=idr[i-1];
    }

    /* For each randomization */
    for (b=1; b<=nsimra; b++) {
	
	/* Compute "used" normal score */
	
	/* for each variable */
	for (i=1; i<=nc; i++) {
	    
	    /* sort the variable */
	    for (j = 1; j <= nla; j++) {
		varR[j-1]=av[j][i];
		indexR[j-1]=j;
	    }
	    
	    /* add a very small value to varR (kind of jitter) */
	    /* first find the smallest lag between successive values */
	    tmp = fabs(varR[1]-varR[0]);
	    
	    for (k = 2; k<=nla; k++) {
		for (l = 1; l<k; l++) {
		    if (fabs(varR[k-1]-varR[l-1]) < tmp) {
			tmp=fabs(varR[k-1]-varR[l-1]);
		    }			    
		}
	    }
	    
	    tmp=tmp/100;
	    
	    for (j = 1; j <= nla; j++) {
		GetRNGstate();
		varR[j-1] = varR[j-1]+ (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
	    }
	    
	    rsort_with_index(varR, indexR, nla);
	    
	    /* normal score for each available point */
	    for (k = 1; k <=nla; k++) {
		nscoav[indexR[k-1]][i] = qnorm((((double) k)/(((double) nla) +1)), 0.0, 1.0, 1, 0);
	    }
	    
	    
	    for (k = 1; k <= nlu; k++) {
		
		/* add a random value to the used points, to avoid ties */
		GetRNGstate();
		res=us[k][i] + (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
		r=1;

		for (j = 1; j <= nla; j++) {
		    
		    /* if the obs is larger than the point, increase r */
		    if (varR[j-1] < res)
			r=j;
		}
		
		tmp3 = (((double) r) + 0.5) / ( ((double) nla) + 1.0);
		/* a qnorm on this, and store the result in nsco */
		nsco[k][i] = qnorm(tmp3, 0.0, 1.0, 1, 0);
		
		
	    }
	    
	}
	
	/*
	  Computes the conditional standardized values, using the method
	  of Ripley (1987, p. 98)
	*/
	/* Just in case there is only one variable */
	if (nc > 1) {
	    for (k = 1; k <= nc; k++) {
		
		l=1;
		for (j=1; j<=nc; j++) {
		    if (j!=k) {
			index[l]=j;
			l++;
		    }
		}
		index[nc]=k;
		
		/* variance covariance matrix of the available points */
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			var[j][l]=0.0;
		    }
		}
		
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			mu1=0.0;
			mu2=0.0;
			
			for (i=1; i<=nla; i++) {
			    mu1 = mu1 + (nscoav[i][index[j]]/((double) nla));
			    mu2 = mu2 + (nscoav[i][index[l]]/((double) nla));
			}
			
			for (i=1; i<=nla; i++) {
			    var[j][l] = var[j][l] + 
				(((nscoav[i][index[j]] - mu1) * 
				  (nscoav[i][index[l]] - mu2))/
				 (((double) nla) - 1.0));
			}
		    }
		}
		
		/* Calculation of Akmo */
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			inv1[j][l] = var[j][l];
		    }
		}
		
		i=invers(inv1, (nc-1), zer, (nc-1));
		
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			Akmo[j][l] = inv1[j][l];
		    }
		}
		
		/* calculation of a */
		for (j = 1; j <= (nc-1); j++)
		    a[j]=var[j][nc];
		
		/* sigkk */
		sigkk=var[nc][nc];
		
		/* for each observation */
		for (i = 1; i <= nlu; i++) {
		    
		    /* Wk */
		    for (j=1; j<= (nc-1); j++)
			Wk[j]=nsco[i][index[j]];
		    
		    /* calculation of the conditionnal average */
		    for (j=1; j<= (nc-1); j++) {
			tmp2[j]=0;
			for (l=1; l<= (nc-1); l++) {
			    tmp2[j] = tmp2[j] + (a[l] *  Akmo[l][j]);
			}
		    }
		    
		    m=0.0;
		    for (l=1; l<= (nc-1); l++) {
			m = m + (tmp2[l] * Wk[l]);
		    }
		    
		    /* Calculation of the conditional variance */
		    s = 0.0;
		    for (l=1; l<= (nc-1); l++) {
			s = s + (tmp2[l] * a[l]);
		    }
		    s = sqrt(sigkk - s);
		    
		    
		    nscob[i][k]=(nsco[i][k] - m)/s;
		}
		
	    }
	} else {
	    nscob[i][1] = nsco[i][1];	    
	}
	

	
	/* ************************
	   Compute the statistics for each individual 
	*/

	/* the number of observations for each animal */
	for (l=1; l<=nid; l++) {
	    e=0;
	    for (i=1; i<=nlu; i++) {
		if (id[i]==l) {
		    e++;
		}
	    }
	    ni[l]=e;
	}
	
	
	/* for each variable */
	for (j=1; j<=nc; j++) {
	    
	    /* for each animal */
	    s=0.0;
	    for (l=1; l<=nid; l++) {
		
		/* get the observations */
		vecalloc(&obs, ni[l]);
		
		k=1;
		for (i=1; i<=nlu; i++) {
		    if (id[i]==l) {
			obs[k]=nscob[i][j];
			k++;
		    }
		}
		
		/* then compute the elements needed for the total variance */
		for (i=1; i<=nsim; i++) {
		    
		    GetRNGstate();
		    tmp= unif_rand();
		    PutRNGstate();
		    for (k = 1; k <= ni[l]; k++) {
			if (tmp < (  ( ((double) k) ) / ((double) ni[l])  )  ) {
			    if (tmp >= (  (((double) k-1)) / ((double) ni[l])  )  ) {
				sammr[i][l]= obs[k];
			    }
			}
		    }
		}
		
		/* and finally, the within class variance */
		m=0.0;
		for (i=1; i<=ni[l]; i++) {
		    m=m+(obs[i] / ((double) ni[l]));
		}
		for (i=1; i<=ni[l]; i++) {
		    s=s+ (((obs[i]-m) * (obs[i]-m)) / 
			  (((double) nlu) - ((double) ni[l])));
		}
		Zi[l]=m;
	    }
	    
	    /* the total variance */
	    vartot=0.0;
	    for (l=1; l<=nsim; l++) {
		/* the mean */
		m=0.0;
		for (i=1; i<=nid; i++) {
		    m=m+ (sammr[l][i]/((double) nid));
		}
		/* the total variance */
		
		for (i=1; i<=nid; i++) {
		vartot=vartot+(((sammr[l][i] - m) * 
				(sammr[l][i] - m))/
			       (((double) nid -1.0) * ((double) nsim)));
		}
	    }
	    sig2=vartot-s;
	    if (sig2<0)
		sig2=0.0;
	    mu[b][((j-1) * 5) + 1]=sammr[1][1];
	    mu[b][((j-1) * 5) + 2]=s;
	    
	    /* value of rho */
	    mu[b][((j-1) * 5) + 3]=sig2/(sig2+s);
	    
	    
	    /* thetai */
	    for (l=1; l<=nid;l++) {
		thetai[l]=sig2 + (s/((double) ni[l]));
	    }
	    
	    /* mu and var(mu) */
	    for (l=1; l<=nid;l++) {
		mu[b][((j-1) * 5) + 4]=mu[b][((j-1) * 5) + 4]+(Zi[l]/thetai[l]);
		mu[b][((j-1) * 5) + 5]=mu[b][((j-1) * 5) + 5]+(1.0/thetai[l]);
	    }
	    mu[b][((j-1) * 5) + 5]=1/mu[b][((j-1) * 5) + 5];
	    mu[b][((j-1) * 5) + 4]=mu[b][((j-1) * 5) + 4] * mu[b][((j-1) * 5) + 5];
	    
	    freevec(obs);
	
	}
	
	/* end of randomization */
    }
    
    k=0;
    for (i=1; i<=nsimra; i++) {
	for (j=1; j<=(nc*5); j++) {
	    resr[k]=mu[i][j];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    usr[k]=nscob[i][j];
	    k++;
	}
    }
 
    /* free memory */

    /* used and available */
    freetab(av);
    freetab(us);
    
    /* the id */
    freeintvec(id);
    
    /* the results */
    freetab(mu);
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	freevec(Wk);
	freetab(Akmo);
	freevec(a);
	freevec(tmp2);
	freeintvec(index);
	freetab(var);
	freetab(inv1);
	freetab(zer);
    }
    freetab(nscob);
    
    /* elements required for the calculations
       of the statistics */
    freetab(sammr);
    freeintvec(ni);
    freevec(Zi);
    freevec(thetai);

    /* elements required for used normal scores */
    freetab(nsco);

    /* elements required for available normal scores */
    freetab(nscoav);


}











/* Extension of the metod for designs I */

void engen2008Ir(double *avr, double *usr, int *nliga, int *nligu, 
		int *ncol, double *resr, int *nsimrar)
{
    /* declaration of variables */
    double **av, **us, **nsco, *varR, tmp, res, **var, mu1, mu2, **Akmo;
    double **inv1, **zer, *a, sigkk, *Wk, *tmp2, m, s, **nscob, **nscoav;
    double tmp3, **mu;
    int nla, nlu, nc, i, j, k, r,l, *index;
    int *indexR, nsimra, b;
    
    /* memory allocation */
    nc=*ncol;
    nla=*nliga;
    nlu=*nligu;
    nsimra=*nsimrar;
    r=1;
    tmp=0.0;
    tmp3=0.0;
    b=1;

    /* used and available */
    taballoc(&av, nla, nc);
    taballoc(&us, nlu, nc);
        
    /* the results */
    taballoc(&mu, nsimra, (nc * 2) );
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	vecalloc(&Wk, (nc-1));
	taballoc(&Akmo, nc-1, nc-1);
	vecalloc(&a, (nc-1));
	vecalloc(&tmp2, (nc-1));
	vecintalloc(&index, nc);
	taballoc(&var, nc, nc);
	taballoc(&inv1, nc-1, nc-1);
	taballoc(&zer, nc-1, nc-1);
    }
    taballoc(&nscob, nlu, nc);


    /* elements required for used normal scores */
    taballoc(&nsco, nlu, nc);
    
    /* elements required for available normal scores */
    taballoc(&nscoav, nla, nc);

    /* elements used in the R API */
    varR = (double *) R_alloc(nla, sizeof(double));
    indexR = (int *) R_alloc(nla, sizeof(int));



        
    /* R -> C objects */
    k=0;
    for (i=1; i<=nla; i++) {
	for (j=1; j<=nc; j++) {
	    av[i][j]=avr[k];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    us[i][j]=usr[k];
	    k++;
	}
    }
    

    /* For each randomization */
    for (b=1; b<=nsimra; b++) {
	
	/* Compute "used" normal score */
	
	/* for each variable */
	for (i=1; i<=nc; i++) {
	    
	    /* sort the variable */
	    for (j = 1; j <= nla; j++) {
		varR[j-1]=av[j][i];
		indexR[j-1]=j;
	    }
	    
	    /* add a very small value to varR (kind of jitter) */
	    /* first find the smallest lag between successive values */
	    tmp = fabs(varR[1]-varR[0]);
	    
	    for (k = 2; k<=nla; k++) {
		for (l = 1; l<k; l++) {
		    if (fabs(varR[k-1]-varR[l-1]) < tmp) {
			tmp=fabs(varR[k-1]-varR[l-1]);
		    }			    
		}
	    }
	    
	    tmp=tmp/100;
	    
	    for (j = 1; j <= nla; j++) {
		GetRNGstate();
		varR[j-1] = varR[j-1]+ (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
	    }
	    
	    rsort_with_index(varR, indexR, nla);
	    
	    /* normal score for each available point */
	    for (k = 1; k <=nla; k++) {
		nscoav[indexR[k-1]][i] = qnorm((((double) k)/(((double) nla) +1)), 0.0, 1.0, 1, 0);
	    }
	    
	    
	    for (k = 1; k <= nlu; k++) {
		
		/* add a random value to the used points, to avoid ties */
		GetRNGstate();
		res=us[k][i] + (((unif_rand() * 2)-1) * tmp);
		PutRNGstate();
		r=1;

		for (j = 1; j <= nla; j++) {
		    
		    /* if the obs is larger than the point, increase r */
		    if (varR[j-1] < res)
			r=j;
		}
		
		tmp3 = (((double) r) + 0.5) / ( ((double) nla) + 1.0);
		/* a qnorm on this, and store the result in nsco */
		nsco[k][i] = qnorm(tmp3, 0.0, 1.0, 1, 0);
		
		
	    }
	    
	}
	
	/*
	  Computes the conditional standardized values, using the method
	  of Ripley (1987, p. 98)
	*/
	/* Just in case there is only one variable */
	if (nc > 1) {
	    
	    /* for each variable */
	    for (k = 1; k <= nc; k++) {
		
		l=1;
		for (j=1; j<=nc; j++) {
		    if (j!=k) {
			index[l]=j;
			l++;
		    }
		}
		index[nc]=k;
		
		/* variance covariance matrix of the available points */
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			var[j][l]=0.0;
		    }
		}
		
		for (j=1; j<=nc; j++) {
		    for (l=1; l<=nc; l++) {
			mu1=0.0;
			mu2=0.0;
			
			for (i=1; i<=nla; i++) {
			    mu1 = mu1 + (nscoav[i][index[j]]/((double) nla));
			    mu2 = mu2 + (nscoav[i][index[l]]/((double) nla));
			}
			
			for (i=1; i<=nla; i++) {
			    var[j][l] = var[j][l] + 
				(((nscoav[i][index[j]] - mu1) * 
				  (nscoav[i][index[l]] - mu2))/
				 (((double) nla) - 1.0));
			}
		    }
		}
		
		/* Calculation of Akmo */
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			inv1[j][l] = var[j][l];
		    }
		}
		
		i=invers(inv1, (nc-1), zer, (nc-1));
		
		for (j = 1; j <= (nc-1); j++) {
		    for (l = 1; l <= (nc-1); l++) {
			Akmo[j][l] = inv1[j][l];
		    }
		}
		
		/* calculation of a */
		for (j = 1; j <= (nc-1); j++)
		    a[j]=var[j][nc];
		
		/* sigkk */
		sigkk=var[nc][nc];
		
		/* for each observation */
		for (i = 1; i <= nlu; i++) {
		    
		    /* Wk */
		    for (j=1; j<= (nc-1); j++)
			Wk[j]=nsco[i][index[j]];
		    
		    /* calculation of the conditionnal average */
		    for (j=1; j<= (nc-1); j++) {
			tmp2[j]=0;
			for (l=1; l<= (nc-1); l++) {
			    tmp2[j] = tmp2[j] + (a[l] *  Akmo[l][j]);
			}
		    }
		    
		    m=0.0;
		    for (l=1; l<= (nc-1); l++) {
			m = m + (tmp2[l] * Wk[l]);
		    }
		    
		    /* Calculation of the conditional variance */
		    s = 0.0;
		    for (l=1; l<= (nc-1); l++) {
			s = s + (tmp2[l] * a[l]);
		    }
		    s = sqrt(sigkk - s);
		    
		    
		    nscob[i][k]=(nsco[i][k] - m)/s;
		}
		
	    }
	} else {
	    /* if one variable, no correction */
	    nscob[i][1] = nsco[i][1];	    
	}
	

	
	/* ************************
	   Compute the statistics
	*/

	/* for each variable */
	for (j=1; j<=nc; j++) {
	    
	    m=0.0;
	    s=0.0;
	    	    
	    /* The preference */
	    for (i=1; i<=nlu; i++) {
		m = m + (nscob[i][j]);
	    }
	    
	    /* the variance of the scores */
	    for (i=1; i<=nlu; i++) {
		s=s+ (((nscob[i][j]-m) * (nscob[i][j]-m)) / 
		      (((double) nlu) - 1));
	    }
	    
	    /* variance of the mean */
	    s=s/((double) nlu);
	    
	    mu[b][((j-1) * 2) + 1]=m;
	    mu[b][((j-1) * 2) + 2]=s;
	}
	
	/* end of randomization */
    }
    
    k=0;
    for (i=1; i<=nsimra; i++) {
	for (j=1; j<=(nc*2); j++) {
	    resr[k]=mu[i][j];
	    k++;
	}
    }
    
    k=0;
    for (i=1; i<=nlu; i++) {
	for (j=1; j<=nc; j++) {
	    usr[k]=nscob[i][j];
	    k++;
	}
    }
 
    /* free memory */

    /* used and available */
    freetab(av);
    freetab(us);
    
    
    /* the results */
    freetab(mu);
    
    /* elements required for the calculations of 
       normalized scores */
    if (nc > 1) {
	freevec(Wk);
	freetab(Akmo);
	freevec(a);
	freevec(tmp2);
	freeintvec(index);
	freetab(var);
	freetab(inv1);
	freetab(zer);
    }
    freetab(nscob);
    

    /* elements required for used normal scores */
    freetab(nsco);

    /* elements required for available normal scores */
    freetab(nscoav);


}


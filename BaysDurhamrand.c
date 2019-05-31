#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*Function prototype */
float ran1(long *idum);
float gasdev(long *idum);
void lowpassnoise(int length, long seed, double out[]);


float gasdev(long *idum)
/*Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
as the source of uniform deviates.*/
/* This function and the random number generator is taken directly from Numerical Recipes in C, 2nd Ed.,
    Cambridge University press, by William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
    ISBN 0-521-43108-5 */
{
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    if (*idum < 0) iset=0; /*Reinitialize. */
        if (iset == 0) { /*We don't have an extra deviate handy, so do */
        do {
            v1=2.0*ran1(idum)-1.0; /* pick two uniform numbers in the square ex*/
            v2=2.0*ran1(idum)-1.0; /* tending from -1 to +1 in each direction, */
            rsq=v1*v1+v2*v2; /* see if they are in the unit circle, */
        } while (rsq >= 1.0 || rsq == 0.0); /*and if they are not, try again.*/
        fac=sqrt(-2.0*log(rsq)/rsq);
        /*Now make the Box-Muller transformation to get two normal deviates. Return one and
        save the other for next time.*/
        gset=v1*fac;
        iset=1; /*Set flag.*/
        return v2*fac;
        } else { /*We have an extra deviate handy,*/
        iset=0; /*so unset the flag, */
        return gset; /*and return it.*/
    }
}



float ran1(long *idum)
/* "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
successive deviates in a sequence. RNMX should approximate the largest floating value that is
less than 1. */

/* This algorithm is good for random arrays less than 100 million */
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0 || !iy) { /*Initialize.*/
if (-(*idum) < 1) *idum=1; /*Be sure to prevent idum = 0. */
else *idum = -(*idum);
for (j=NTAB+7;j>=0;j--) { /* Load the shue table (after 8 warm-ups). */
k=(*idum)/IQ;
*idum=IA*(*idum-k*IQ)-IR*k;
if (*idum < 0) *idum += IM;
if (j < NTAB) iv[j] = *idum;
}
iy=iv[0];
}
k=(*idum)/IQ; /* Start here when not initializing. */
*idum=IA*(*idum-k*IQ)-IR*k; /* Compute idum=(IA*idum) % IM without overflows */
if(*idum < 0) *idum += IM; /* by Schrage's method. */
j=iy/NDIV; /* Will be in the range 0..NTAB-1. */
iy=iv[j]; /* Output previously stored value and refll the */
iv[j] = *idum; /* shue table. */
if ((temp=AM*iy) > RNMX) return RNMX; /* Because users don't expect endpoint values.*/ 
else return temp;
}

void lowpassnoise(int length, long seed, double out[])
{
    /*This function generates normally distributed random variables and then low-pass filters it
    with a 4th order Butterworth filter */
    int n;
    double *x,
    b[]={0.00011138, 0.00044552, 0.00066829, 0.00044552, 0.00011138}, a[]={1, -3.4259, 4.4368, -2.5712, 0.5621};
    /*Parameters for 4th order butterworth filter with 500Hz cutoff and sampling frequency = Fs/2, where Fs=1/delt/tau_m=1/0.01/0.007*/
    /* that leads to a ratio of cutoff to nyquist of 0.07 */
    /*IF YOU CHANGE TAU_M OR DELT, YOU NEED TO LOOK UP DIFFERENT PARAMETERS FOR BUTTERWORTH FILTER (FOR EX. MATLAB HAS THEM)*/

    
    /*Initializing vector to place normally distributed numbers in*/
    x=(double *) calloc(length,sizeof(double)); 
    
    for(n=0;n<length;n++)
    {
        x[n]=gasdev(&seed); /*Gets Gaussian random numbers*/
        if(n<4) out[n]=x[n];    
        else out[n]=b[0]*x[n] + b[1]*x[n-1]+b[2]*x[n-2]+b[3]*x[n-3]+b[4]*x[n-4] - a[1]*out[n-1] - a[2]*out[n-2] - a[3]*out[n-3] -a[4]*out[n-4];
    }           /*  -------------Feedforward filter parameters----------------    --------------Feedback filter parameters------------------ */               
    free(x); 
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LIFDAPC.c"
#include "BaysDurhamrand.c"
#define pi 3.141593

int
main(void)
{
    /*---Local variables---*/
    long seed;
    int i,j, GL, numw, nbins, ISIbins;
    double *electroRinput, *lowpassnoisevector, *ISIedges, *PSTHedges, ISIbinwidth,Fs, ISI, *weightsbefore, *weightsafter, numberofT, spmodT;
    
    /* OUTPUTS*/
    double *sptime, *rec, avgw=0, *PSTHhist, *PSTHphase, *ISIhist, firing, burst2,burst4; /* Dynamically allocated output arrays*/
    int spsize, nspikes, n2bursts, n4bursts;
    FILE *outpsth, *outparam, *outISI;
    
    /*---Simulation Parameters */
    double delt=0.01, endtime=1750, freq; /* in seconds */
    int maxj;

    /*---Feedforward parameters---*/
    double tau_m= 0.007, fcut1=500, sigma=0.759, I=0.759*0.759, kappa=0.39;

     /*---Feedback parameters---*/
    double eta=0.0036, eta2, Lambda;
    /*if 2- and 4-spike burst, optimal fit is*/
    double wmax=1.5, tau_w=980, g;
    /*if 2-spike burst only, optimal fit is below, with g=0.87*Lambda; */
    /* double wmax=1, tau_w= 4900,g; */
    
    

    
        
   
    /*---Ask user for frequency---*/
    printf("what is the AM frequency?\n");
    scanf("%lf",&freq);
    getchar();
    
    printf("if global, enter 1\n");
    scanf("%d",&GL);
    getchar();

    /*get seed*/
    printf("enter a seed (negative integer) \n");
    scanf("%d",&seed);
    if(seed>0)
    seed=-1*seed;
    
    getchar();

    numw= floor(1/freq/0.0025);

   

    /*electroreceptor adaptation*/
    if(freq==0.5)
        kappa=0.25;
    if(freq==1)
        kappa=0.27;
    if(freq==2)
        kappa=0.31;

   
   

  


    weightsbefore =(double *) calloc(numw,sizeof(double)); 
    weightsafter =(double *) calloc(numw,sizeof(double)); 
    for(i=0;i<numw;i++)
    {
        weightsbefore[i]=wmax;
    }


    if (GL==1)
    {
        Lambda = 1;
        maxj=3;
    }
    else
    {
        Lambda=0;
        maxj=1;
    }

    g=1.44*Lambda;
    
    /*eta decreasing at high freqs to simulate granule cell burst decay */
/*
    if (freq ==12)
        eta = eta/2;
    if (freq ==16)
        eta=eta/4;
    if (freq==20)
        eta=0.0003;
    if (freq==32)
        eta=0.00015;

    eta2=eta/2;
    */

    /*Lambda decreasing at high freqs model*/
/*   if (freq ==12)  Lambda =0.5; 
    if (freq ==16)  Lambda =0.27;
    if (freq==20)   Lambda =0.05;
    if (freq==32)   Lambda =0;
*/

/*---Initialization*/

    spsize= endtime/delt/tau_m/100;
    sptime=             (double *) calloc(spsize, sizeof(double));
    rec=                (double *) calloc(spsize, sizeof(double));
    electroRinput =     (double *) calloc(endtime/delt/tau_m, sizeof(double));
    lowpassnoisevector= (double *) calloc(endtime/delt/tau_m, sizeof(double));

    Fs=1/delt/tau_m; /* sampling rate in real time*/ 

    
    /*---Simulation loop---*/
    for(j=0;j<maxj;j++)
    { 
        /*This function generates lowpass filtered gaussian noise with cutoff frequency 500Hz, assuming Fs =1/0.01/0.007
        If you change tau_m or delt, you NEED to change the lowpassnoise function*/
        lowpassnoise(endtime/delt/tau_m,seed,lowpassnoisevector);
        
        /*feedforward input*/
        for(i=0;i<endtime/delt/tau_m; i++)
        {
            electroRinput[i]=I+kappa*sin(-2*pi*freq*i*delt*tau_m)+sigma/sqrt(fcut1*2/Fs)*lowpassnoisevector[i];
            if(electroRinput[i]<0)electroRinput[i]=0; /*rectification*/
        }
    
        LIFDAPmodel(electroRinput,endtime/delt/tau_m,freq,weightsbefore,numw,g,Lambda,
        eta,eta2,tau_m,tau_w,wmax,sptime,rec,weightsafter, &nspikes, &n2bursts, &n4bursts); 
   
        for(i=0;i<numw;i++) weightsbefore[i]=weightsafter[i];
        
       

    }
    
    /*Printing average firing and burst rates. If global, also printing average weights*/
    burst4=n4bursts/endtime; /*4-sp burst rate in seconds*/
    burst2=n2bursts/endtime; /*2-sp burst rate in seconds*/
    firing=nspikes/endtime;

    printf("average spike rate is %f\n",firing);
    printf("average 2-spike burst rate is %f\n",burst2);
    printf("average 4-spike burst rate is %f\n",burst4);

    if(GL==1)
    {
        for(i=0;i<numw;i++) avgw += weightsafter[i]/numw;
        printf("The average weight value was %f\n",avgw);
    }
    getchar();

    outparam=fopen("params.txt","w");
        fprintf(outparam,"average spike rate is %f\n",firing);
        fprintf(outparam,"average 2-spike burst rate is %f\n",burst2);
        fprintf(outparam,"average 4-spike burst rate is %f\n",burst4);
        fprintf(outparam,"The average weight value was %f\n",avgw);
     
     fclose(outparam);

 

    
    /*Histogramming firing times. nbins is the number of bins period is divided into*/
    nbins=20; /*number of bins of PSTH over period*/
    ISIbins=50; /*number of ISI bins*/
    ISIbinwidth=4.0; /*Width of each ISI histogram bin in milliseconds*/

    numberofT=endtime*freq; /*number of periods simulated*/
    PSTHedges=(double *) calloc(nbins+1, sizeof(double));
    ISIedges= (double *) calloc(ISIbins+1, sizeof(double));
    PSTHhist=(double *) calloc(nbins, sizeof(double));
    PSTHphase=(double *) calloc(nbins, sizeof(double));
    ISIhist= (double *) calloc(ISIbins, sizeof(double));

    for(i=0; i<nbins;i++)
    {
        PSTHhist[i]=0.0; /*Initializing histogram */
        PSTHedges[i]=i/freq/nbins; /*Dividing period into equal bins for histogramming*/
        PSTHphase[i]=pi/nbins +2*pi/nbins*i;
    }
    PSTHedges[nbins]=1/freq; /*adding last bin*/

    for(i=0;i<ISIbins;i++)
    {
        ISIhist[i]=0;
        ISIedges[i]=ISIbinwidth*i; /*Each bin covers 4 ms */
    }
    ISIedges[ISIbins+1]=ISIbins*ISIbinwidth; /*longest ISI shown in histogram. (200ms with these values)*/

    for(i=0;i<nspikes;i++)
    {
        spmodT=fmod(sptime[i+4]*tau_m,1/freq); /*finding firing time modulo the period (firing times were in units of tau_m)*/
        /*note that sptimes is offset by 4 in LIFDAP so that 4-spike burst searching can start immediately*/

        for(j=0;j<nbins;j++)
        {
            if(PSTHedges[j]<=spmodT && spmodT<PSTHedges[j+1]) PSTHhist[j]+=1.0/numberofT*freq*nbins; /*binning spike times and converting to Hz*/
        } /*when hist has totaled all spikes then multiplying by 1/numberofT will equal average number of spikes per bin*/
            /*To convert to average firing rate in bin, divide by duration of bin = 1/freq/nbins     */

        if(i<nspikes-1)
        {
            ISI=(sptime[i+4+1]-sptime[i+4])*tau_m*1000; /*ISIs are in milliseconds now!*/
            for(j=0;j<ISIbins;j++)
            {
                if(ISIedges[j]<=ISI && ISI<ISIedges[j+1]) ISIhist[j]+=1.0/nspikes/ISIbinwidth; /*Normalized to 1*/ 
            }

        }
    }

    /*print PSTH to file*/
    outpsth=fopen("PSTH.txt","w");
    for(i=0;i<nbins;i++)
    
     fprintf(outpsth, "%f \t %f \n",PSTHphase[i],PSTHhist[i]); /*So PSTHs are always plotted against phase.
                                                                 That way, different frequencies can be easily plotted against each other.*/
    
     fclose(outpsth);

   /*Print ISI to file*/
    outISI=fopen("ISI.txt","w");
    for(i=0;i<ISIbins;i++)
    
     fprintf(outISI, "%f \t %f \n",ISIedges[i]+ISIbinwidth/2,ISIhist[i]);
    
     fclose(outISI);
   
return(0);

}

/*

%histogramming the firing times
npts=20; %20 bins in histogram
sptimes=A(A>0)*tau_m; %all the spike times (A is in units of tau_m)
spper= rem(sptimes,abs(1/freq));  %finding spike times mod period of AM

edges=linspace(0,abs(1/freq),npts+1); 
n=histc(spper,edges);
numper=endtime*freq; %number of periods avgd over to get mean PSTH
n=n/numper*npts*freq; %n/numper = avg #  spikes in the bin in 1 period; 
%T/npts =1/(f*npts) = length of each bin 



%ISI stuff
ISI=sptimes(2:end)-sptimes(1:end-1);
ISI=ISI*1000; %ISI now in ms
ISIn=histc(nearest(ISI*10)/10,0:200); %was 501
ISIn=ISIn/sum(ISIn(1:end-1));


*/



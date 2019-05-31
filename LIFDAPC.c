#include <math.h>


#if !defined(MAX)
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#endif

/*Function Prototypes*/
void LIFDAPmodel(double signal[], int length_signal, double f, double weight[],
       int length_weight, double g, double Lambda, double eta, double eta2, double tau_m, double tau_w, double wmax, double sptime[],
       double rec[], double weightafter[], int *nspikes, int *n2bursts, int *n4bursts);                   

void governingloop(double sptime[], double bursttime4[], double bursttime2[],
        double w[], double signal[], double weight[], double rec[],double f, double g, double Lambda, 
        double eta, double eta2, int numw, int m, double tau_m, double tau_w, double wmax,int *nspikes, int *n2bursts,int *n4bursts);

/*This calling function only exists to call the time loop subroutine. This parallels the MATLAB calling function in LIFDAP.c */


void LIFDAPmodel(double signal[], int length_signal, double f, double weightbefore[], int length_weight, double g, double Lambda,
double eta, double eta2, double tau_m, double tau_w, double wmax, double sptime[], double rec[], double weightafter[], int *nspikes, int *n2bursts, int *n4bursts)                   
{

  double  *bursttime2, *bursttime4;
  int spsize; 

    spsize=length_signal/100;

    bursttime2= (double *) calloc(spsize, sizeof(double));
    bursttime4= (double *) calloc(spsize, sizeof(double));
    
    /* Call the LIF-DAP timeloop subroutine. */
    governingloop(sptime,bursttime4,bursttime2, weightafter,signal,weightbefore,rec,f,g,Lambda,eta,eta2,
     length_weight,length_signal,tau_m, tau_w, wmax, nspikes, n2bursts, n4bursts);
    /*spttime, bursttime4, bursttime2, w and rec are outputs that are arrays
    signal and weight are inputs that are arrays
    the rest are individual inputs */


}


/*-----------------------------------------------------------------------------------------------*/

void governingloop(double sptime[], double bursttime4[], double bursttime2[],
        double w[], double signal[], double weight[], double rec[],double f, double g, double Lambda, 
        double eta, double eta2, int numw, int m, double tau_m, double tau_w, double wmax, int *nspikes, int *n2bursts,int *n4bursts)
{
    /*Parameters from Noonan et al. 2003*/
 double A=0.15*4, B=2, alpha=20, 
        beta=0.35, D=0.1, E=3.5, tauref=0.1,taudend=50., b=0, 
        somawidth=0.05*4, dendwidth=1.0, taudecay=1., thresh=1.0;
     
 double delt = 0.01, realt, avgw, Lwidth, Lwidth2, Bdef4,Bdef2, nper, burstT, *pfspike, *L, v=0.025, tref=0, Dxh, Dwh, Dyh, Dsh, Dx=0, Dy=0, Ds=0, Dw=0; 
 int i, k, j, reci=0, index=4, index4=0, index2=0, count4=4, count2=2,  countr=0;
    

  
  /*PF initialization */
    pfspike= (double *) calloc(3*numw, sizeof(double));
    L=(double *) calloc(3*numw, sizeof(double));
   
 /*Burst definition and width*/   
    Bdef2=0.015/tau_m; /*defining a 2-spike burst (time is in units of tau_m within code so 15 ms = 0.015/tau_m)*/      
    Bdef4=3.0*Bdef2; /*4-spike burst has 3 ISIs so 3 x 2-spike burst definition */
    Lwidth=0.1/tau_m; /*100 ms 4-spike burst learning rule width, from experimental data*/
    Lwidth2=0.01/tau_m; /*10 ms 2-spike STDP width*/
    tau_w=tau_w/tau_m;
    f=f*tau_m; 
    nper=1/f/delt; /*number of timesteps in 1 period*/

   /*Mapping given weight distribution to output weight matrix w */
   for(i=0;i<numw;i++){
       w[i]=weight[i];
   }
   
   /* pfspike gives start times (i.e. firing times) of each PF over 3 periods  */
   for (i=0;i<3*numw;i++){
       pfspike[i]=ceil(nper/numw*i)*delt; /* nper/numw*delt= T/numw = time span of each segment (i.e. how long is each weight "on" for = 2.5ms) */
   }
    
   
    /* Initializations*/
    sptime[0]=-100; sptime[1]=-100; sptime[2]=-100;    sptime[3]=-100; 
    
   
    
    /*----Time loop----*/
    for (i=0;i<m;i++){
        realt=delt*i; /*realt in units of tau_m, so not real time per say*/
        k= (int) floor(fmod(numw*realt/nper/delt,numw));
        /*k is an integer that increases stepwise and signals the start of a new PF */
        /*so for 0-2.5ms, k=0, for 2.5-5ms, k=1, etc., and loops back to k=0 when period is over (hence modf of numw)*/
        
        
        
        /*if realt is greater than the refractory period*/     
        if(realt>tref){
            /*GOVERNING EQUATION*/           
            v= v+delt*(-v + signal[i] + Lambda*(w[k]-g*v));
            /*k is used to change the index of the weight that is active at a given time*/
            /*Note absence of DAP*/
                
                
        
            
            
            /* Do this if ISI beyond dendritic ref. period*/
            if(taudend<(sptime[index-1]-sptime[index-2])){ 
                if (dendwidth*Dx-somawidth*Ds > 0){  /*DAP is rectified*/  
                    v=v+delt*alpha*(dendwidth*Dx-somawidth*Ds); /*DAP is added*/
                }
            }
        
        }/*end of if realt > tref */
        
        /*if V >threshold, a spike is fired and a burst might be recorded*/
        if(v>thresh) {
                v=0; /*v reset to 0*/
                tref=realt+tauref+delt/2; /*refractory period is updated*/
                sptime[index]=realt; /*spike is recorded*/
                index++; /*index now moves to vacant position*/
                count4--; /*so each 4-spike burst has 4 unique spikes*/
                count2--; /* so each 2-spike burst has 2 unique spikes */
                
                /*DAP parameters are updated*/
                b=b+A+ B*b*b; /*updating b*/
                taudend=D+E*b; /*updating dendritic refractory period*/
                dendwidth=beta*b;
                Dy=Dy+1/dendwidth/dendwidth;
                Dw=Dw+1/somawidth/somawidth;
   
                /*if the last spike occurred within Bdef2 of this spike and count2<1, then record a 2-sp burst*/
                if((realt-sptime[index-2]<Bdef2)&& (count2<1)) { /*count2 makes sure that bursts don't share spikes */
                    bursttime2[index2]=sptime[index-2];/* tracks the 1st spike in burst (hence "index-2")*/
                    index2++; /*2sp burst index moves up 1*/
                    count2=2; /*count reset*/

                    /*learning 2sp rule*/
                    burstT= fmod(sptime[index-2],nper*delt)+nper*delt; 
                    /*burstT = time of SP burst, mod the period of AM (i.e. at 4 Hz, 550ms = 50 ms after start of period) + 1 period*/
                    /*to make sure that a burst at the end of a period affects weights at the beginning of the next cycle,
                      and same with a burst at the beginning of a period affecting weights at the end of the last cycle,
                      pfspike has PF "firing" times for 3 periods and burstT adds a period (i.e. +nper*delt) to the burst time*/
                  
                    /*Also, since I know exactly when PFs will fire in the future, 
                     I apply the learning rule both pre-post and post-pre when the SP cell fires*/
                    
                    
                    for (j=0;j<3*numw;j++){
                        L[j]=1-pow((pfspike[j]-burstT)/Lwidth2,2); /*Quadratic Learning rule for each PF time */
                        if(L[j]<0){L[j]=0;} /*rectification of the learning rule (so it's strictly inhibitory) */
                    }
                    for (j=0;j<numw;j++){
                        w[j]=w[j]-eta2*w[j]*(L[j]+L[j+numw]+L[j+2*numw]); /*weights updated*/
                        if(w[j]<0){w[j]=0;} /*depression at that weight's segment from each of the 3 periods looked at is added together*/
                      
                    } /*for 2-spike bursts, Lwidth2 is small, so L <0 --> L=0 often*/

                }/*End of if 2-spike burst occurred*/
                
                
             /*if the time between this spike and the 4th last spike is less than Bdef4, record a 4-spike burst*/   
             if((realt-sptime[index-4]<Bdef4)&& (count4<1)) { 
                    bursttime4[index4]=sptime[index-4];/* tracks the 1st spike in burst*/
                    index4++;
                    count4=4; /*no overlapping 4 sp bursts*/
                    count2=2; /* so 2sp burst can't use last spike in 4 sp burst*/
                 
                   
                    /*since weights change immediately, once a 4-spike burst is identified,
                      a 2-spike burst has likely just occurred and must be removed
                      (so that a 4-spike burst is not mistakenly double counted as also having 2-spike bursts in it) */
                    
                     /*UNLEARNING LOOP: 2sp bursts within the 4sp burst*/
                    while(bursttime4[index4-1]-bursttime2[index2-1] < delt){
                     
                        burstT= fmod(bursttime2[index2-1],nper*delt)+nper*delt;
                        for (j=0;j<3*numw;j++){
                            L[j]=1-pow((pfspike[j]-burstT)/Lwidth2,2);
                            if(L[j]<0){L[j]=0;}
                        }
                        for (j=0;j<numw;j++){
                            w[j]=w[j]/(1-eta2*(L[j]+L[j+numw]+L[j+2*numw]));
                            if(w[j]<0){w[j]=0;}
                        }
                        /*this finds the effect of 2-sp burst that happened and does the inverse operation
                          Technically, the weights have changed since the burst because of potentiation rule, 
                         but it is negligible (time elapsed ~45 ms << tau_w = 980s) */
                        index2--; /*record of 2sp burst erased*/
                     
                        countr++;
                    } /* repeat unlearning loop until no 2sp burst in last 4 spikes (i.e. could be 0, 1, or 2 bursts)*/
                    
                    /*UNLEARNING LOOP: 2sp burst that used the 1st spike in 4sp burst as its last spike*/
                     if (bursttime2[index2-1]==sptime[index-5]){
                        burstT= fmod(bursttime2[index2-1],nper*delt)+nper*delt;
                        for (j=0;j<3*numw;j++){
                            L[j]=1-pow((pfspike[j]-burstT)/Lwidth2,2);
                            if(L[j]<0){L[j]=0;}
                        }
                        for (j=0;j<numw;j++){
                            w[j]=w[j]/(1-eta2*(L[j]+L[j+numw]+L[j+2*numw]));
                            if(w[j]<0){w[j]=0;}
                        }
                      
                        index2--;
                       
                        countr++;
                    } /* 5th spike unlearning loop*/

                   
                    /*learning 4sp rule*/
                    burstT= fmod(sptime[index-4],nper*delt)+nper*delt; /*range= [T,2T) */
                   
                    /*with Lwidth4 being large compared to T at high AM freqs, using 3 periods means that
                     sometimes one PF will be affected multiple times by 1 burst (i.e. if a PF bursts 50ms 
                     before and 50ms after a SP cell bursts, then the PF will be depressed by the sum of both). 
                     This effect is limited to 3 periods, so at 20 Hz and esp. at 32 Hz, the effect of 4-spike 
                     bursts are clipped at the ends. This was for computational simplicity, but it is unknown 
                     how this situation is resolved in vivo anyway. */
                    
                    for (j=0;j<3*numw;j++){
                        L[j]=1-pow((pfspike[j]-burstT)/Lwidth,2); /*Learning rule*/
                        if(L[j]<0){L[j]=0;}
                    }
                    for (j=0;j<numw;j++){
                        w[j]=w[j]-eta*w[j]*(L[j]+L[j+numw]+L[j+2*numw]); /*weight update*/
                        if(w[j]<0){w[j]=0;}
                    }
                    
                    
             } /*end of if 4sp burst... */
           
             
              
               
                
        } /*end of if fired...*/
            
        /*Dendritic alpha f'n: how to code it dynamically*/
        Dxh=Dx+delt*Dy; /*%D for DAP = dendritic after-polarization*/
        Dyh=Dy+delt*(-Dx/(dendwidth*dendwidth)-2*Dy/dendwidth); 
        Dx=Dxh;
        Dy=Dyh;

        /*Somatic alpha f'n */
        Dsh=Ds+delt*Dw; 
        Dwh=Dw+delt*(-Ds/(somawidth*somawidth)-2*Dw/somawidth); 
        Ds=Dsh;
        Dw=Dwh;
        
        b=b +delt*(-b/taudecay); /*b dynamically decays*/    
        
        
        
        
        /*Potentiation rule*/
        for (j=0;j<numw;j++){
            w[j]=w[j]+delt/tau_w*(wmax-w[j]);
        }
      
    }/*----end of time loop----*/

    *nspikes  = index-4; /*index initializes at 4 so 4spike burst search can start immediately*/
    *n2bursts = index2;
    *n4bursts = index4;
    
} /*end of governingloop subfunction*/




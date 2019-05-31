
function[firing burst2 burst4]=CLSglobal(freq,GL)



%selects local (GL= not one) and global (GL==1)
%under global, default is no correction at high frequencies;
% to simulate granule cell burst rate decay using either eta or lambda
% decay, uncomment the appropriate section

if GL==1
    Lambda = 1; maxj=3;
else
    Lambda=0; maxj=1;
end



  

%---Simulation Parameters

delt=0.01;% 
endtime=1750; %end time in seconds;

numw=floor(1/freq/0.0025);

%---Feedforward parameters---
tau_m= 0.007; %membrane time constant in seconds
fcut1=500;
sigma=0.759;
I=0.576;
kappa=0.39;%0.39;

%electroreceptor adaptation
if freq==0.5
    kappa=0.25;
end
if freq==1
    kappa=0.27;
end
if freq==2
    kappa=0.31;
end

%---Feedback parameters---
%if 2- and 4-spike burst, optimal fit is
wmax=1.5;  tau_w =980; g=1.44*Lambda;% 1.66 if 2-spike only
%if 2-spike burst only, optimal fit is
%wmax=1.; tau_w= 4900; g=0.87*Lambda; %490, 0.95 for single spikes
weights=wmax*ones(1,numw);



eta=0.0036;
%{
%eta decreasing at high freqs model
if freq ==12
    %eta = eta/2;
    eta=eta/4; %2-spike bursts only 
    %eta=eta/3; %for single spikes
end
if freq ==16
    %eta=eta/4;
    eta=eta/10; % 2-spike burst only
    %eta=0.0005; %single spikes 
end
if freq==20
   % eta=0.0003;
    eta=0.00007;    %2-spikes
   %eta=0.0002; %single spikes
end
if freq==32
    %eta=0.00015;
    eta=0.00006; %2-spikes
    %eta=0.00012; %single spikes
end
%}
eta2=eta/2;

%{
%Lambda decreasing at high freqs model
if freq==12
    Lambda =0.5;
end
if freq==16
    Lambda=0.27;
end
if freq==20
    Lambda=0.05;
end
if freq==32
    Lambda=0;
end

%}

%---Initialization
t=(delt*tau_m:delt*tau_m:endtime); %t is in real time
fr=zeros(1,maxj);
Fs=1/delt/tau_m; % sampling rate in real time 
[bf,af]=butter(4, fcut1*2/Fs);
    

%---Simulation loop---
for j=1:maxj
    %---Noise Generation
    gaus=randn(1,endtime/delt/tau_m); 
    lowpass=filter(bf,af,gaus)/sqrt((fcut1*2/Fs));
 
    %feedforward input
     electroRinput=I+kappa*sin(-2*pi*freq*t) +lowpass*sigma;
     electroRinput(electroRinput<0)=0; %rectification
    
    
  
      
    [A,B,C,D,E]=LIFDAPmatlab(electroRinput,freq,weights,g,Lambda,eta,eta2,tau_m,tau_w,delt,wmax); %change  to kappa
   
    
    weights=D;
    fr(j)=length(A(A>0))/endtime;
    
%DO it again maxj times, so transients decay and get data from last epoch
% use maxj=1 for local and maxj=3 for global
end

burst4=length(C(C>0))/endtime; %4-sp burst rate in seconds
burst2=length(B(B>0))/endtime; %2-sp burst rate in seconds
firing=fr(maxj);
   

avgw=mean(D);

%histogramming the firing times
npts=20; %20 bins in histogram
sptimes=A(A>0)*tau_m; %all the spike times (A is in units of tau_m)
spper= rem(sptimes,abs(1/freq));  %finding spike times mod period of AM

edges=linspace(0,abs(1/freq),npts+1); 
n=histc(spper,edges);
numper=endtime*freq; %number of periods avgd over to get mean PSTH
n=n/numper*npts*freq; %n/numper = avg #  spikes in the bin in 1 period; 
%T/npts =1/(f*npts) = length of each bin 

output=n(1:20);

%ISI stuff
ISI=sptimes(2:end)-sptimes(1:end-1);
ISI=ISI*1000; %ISI now in ms
ISIn=histc(nearest(ISI*10)/10,0:200); %was 501
ISIn=ISIn/sum(ISIn(1:end-1));




%
%Plotting
figure
%the PSTH
subplot(2,2,1); plot(pi/npts:2*pi/npts:2*pi-pi/npts,n(1:npts));
hold all
xlabel('PSTH')
xlim([0 6.3])
ylabel('Firing rate (Hz)')
title(['mean f.r. = ', num2str(fr(maxj),'%4.2f'), ' and freq = ',int2str(freq)])

%The weight distribution
subplot(2,2,2);plot(weights)%plot(pi/numw:2*pi/numw:2*pi-pi/numw,weights);
%xlim([0 6.3]);
ylim([0 1.55]);
xlabel('weights')
title(['mean w = ' , num2str(mean(weights),'%5.3f'),' and pk2pk = ', num2str(max(weights)-min(weights),'%4.2f')])

%The ISI histogram
subplot(2,2,3);plot(0.5:199.5,ISIn(1:end-1))
hold all
xlabel('ISI (ms)')
xlim([0 200]);
title(['4-sp bursts = ',num2str(burst4,'%5.3f'),' and 2-sp bursts = ',num2str(burst2,'%5.3f')])

%}




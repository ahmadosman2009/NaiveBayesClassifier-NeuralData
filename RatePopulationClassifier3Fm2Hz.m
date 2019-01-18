  % function [Results] = RatePopulationClassifier3Fm2Hz(Classifier,region,NB,Fm,units)
%   
%       DESCRIPTION     : This function uses a bernoulli bayesian calssifier to estimate the percent correct of single 
%                         trails compared to PSTH templates using spike rate. 
%       in              : All rasters from multiple neurons. Rastergram Data Structure
%                         spet: spike event time 
%                         Fs: sampling rate
%       region          : A1=1 VAF=2 SRAF=3
%       FsHist          : Desired Sampling Rate for cycle PSTH (Hz)
%       NB              : Number of bootstraps to be used. 
%       units           : Number of units to be used in classifier. 
%
%Returned Values
%       Results         : Results from comparing single trials to the PSTH
%                         templates. (input sound,outputsound,bootstrap
%                         number) 
% (C) Ahmad Osman & Monty Escabi, Feb 2017

function [Results] = RatePopulationClassifier3Fm2Hz(Classifier,region,NB,Fm,units)

NU=length(Classifier(region).RastersFm2);  %Total number of units for specified region
N=length(Classifier(region).RastersFm2(1).Units);%Number of sounds
L=length(Classifier(region).RastersFm2(1).Units(1).Trials);%Number of trials 
Fs=Classifier(1).RastersFm2(1).Units(1).Trials(1).Fs;   %spike train sampling rate
Results=zeros(N,N,NB);  %Results matrix containig 0 (miss) or 1 (hit)
Tp=1/Fm;    %Sound Period

for l=1:NB %Bootstraps
    %Define units used for each bootstrap
    Units(l,:)=randsample(NU,units);  
    %Generating Model
    for n=1:N   %soundsl
            %Find model and validation trials            
            for w=Units(l,:)
                %Genearging raster for model trials only
                imod(n,:,w)=randsample(L,L-1)';   %randomly remove 1 trial - model trials
                ival(n,w)=sum(1:L)-sum(imod(n,:,w)); %validation trial
                [RAS]=rasterexpand(Classifier(region).RastersFm2(w).Units(n).Trials(imod(n,:,w)),Fm,4);
                lambda(w,n)= mean(mean(RAS(:,1:8)));    %Average firing rate for model trials - this is our model (Poisson)         
            end

    end
    %Classifying
    for m=1:N   %input sounds
%             LogLike=zeros(1,N);
            for q=1:N   %ouput sounds (model)
                 %Validating for current bootstrap itteration
                 LogLike(q)=0;
                 for w=Units(l,:)
                    [RASv]=rasterexpand(Classifier(region).RastersFm2(w).Units(m).Trials(ival(q,w)),Fm,4);       %Validatoin response - single trial with multiple cycles
                    SpikeCountV=RASv(:,1:8)*Tp;   
                    %Spike count for each validation cycle
                    P=poisspdf(SpikeCountV,lambda(w,q)*Tp);
                    LogLike(q)= LogLike(q)+sum(log10(P)); 
                 end
            end
            
    %Makind decision to maximize log-like
    [hit,index]=max(LogLike);
    Results(m,index,l)=1;
    end  
end 


                 
                 
%             GModel=0;   
%             count = 0;
%             err_count = 0;
%             while count == err_count
%                 try
%                      GModel=gmdistribution.fit(ModelData',3);
%                     if GModel == 0
%                         display('not working')
%                     end
%                 catch display('failed')
%                     err_count = err_count + 1;  
%                 end
%                 count = count + 1;
%             end

            

%                  dx=XX(2)-XX(1)
%                  figure
%                  bar(XX,N/sum(N)/dx)
%                  hold on
%                  plot(XX',pdf(GModel,XX'),'r')
%                                  
%                  X2=randn(1,1024)*.5+3
%                  X1=randn(1,1024);
%                  X=[X1 X2(1:512)]
%                  [N,XX]=hist(X,20)
%                  dx=XX(2)-XX(1)
%                  bar(XX,N/sum(N)*dx)
%                  hold on
%                  plot(XX',pdf(GMModel,XX'),'r')
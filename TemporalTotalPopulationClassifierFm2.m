  % function [Results]=TemporalTotalPopulationClassifier(in,NB,units)
%   
%       DESCRIPTION     : This function uses a bernoulli bayesian calssifier to estimate the percent correct of single 
%                         trails compared to PSTH templates. 
%       in              : All rasters from multiple neurons. Rastergram Data Structure
%                         spet: spike event time 
%                         Fs: sampling rate
%       FsHist          : Desired Sampling Rate for cycle PSTH (Hz)
%       NB              : Number of bootstraps to be used. 
%       units           : Number of units to be used in classifier. 
%
%Returned Values
%       Results         : Results from comparing single trials to the PSTH
%                         templates. (input sound,outputsound,bootstrap
%                         number) 
% (C) Ahmad Osman & Monty Escabi, Feb 2017

% for g=1:50
% in(1).TotalRasters(g).Units=in(1).Rasters(g).Units
% end 
% 
% for g=1:50
% in(1).TotalRasters(g+50).Units=in(2).Rasters(g).Units
% end
% 
% for g=1:50
% in(1).TotalRasters(g+100).Units=in(3).Rasters(g).Units
% end



function [Results]=TemporalTotalPopulationClassifierFm2(Classifier,NB,units)

%units=5; 
FsHist=200;
%W=length(in(region).Rasters); %number of units 
% N=length(Classifier(1).TotalRastersFm2(1).Units);%Number of sounds
% L=length(Classifier(1).TotalRastersFm2(1).Units(1).Trials);%Number of trials 
% Fs=Classifier(1).TotalRastersFm2(1).Units(1).Trials(1).Fs;   %spike train sampling rate
N=length(Classifier(1).TotalRastersFm2fifty(1).Units);%Number of sounds
L=length(Classifier(1).TotalRastersFm2fifty(1).Units(1).Trials);%Number of trials 
Fs=Classifier(1).TotalRastersFm2fifty(1).Units(1).Trials(1).Fs;   %spike train sampling rate


for l=1:NB %Bootstraps
    
for w=1:units
    for n=1:N   %output sounds
                %Find model and validation trials
                imod(n,:,w)=randsample(L,L-1)';   %randomly remove 1 trial - model trials
                ival(n,w)=sum(1:L)-sum(imod(n,:,w)); %validation trial

                %Generate Model Cycle Raster
                %[RASTER]=raster2cycleraster(Classifier(1).TotalRastersFm2(w).Units(n).Trials(imod(n,:,w)),4,1,0,0); % creating cycle raster 
                [RASTER]=raster2cycleraster(Classifier(1).TotalRastersFm2fifty(w).Units(n).Trials(imod(n,:,w)),4,1,0,0); % creating cycle raster 
                [PSTHi(n,:,w)]=cyclepsthmodel(RASTER,FsHist,1000,16); %creating PSTH template
                [PSTHi(n,:,w)]=PSTHi(n,:,w)/sum(PSTHi(n,:,w));
    end
end


for m=1:N   %input sounds
           % LogLike=zeros(1,N);
            for n=1:N   %ouput sounds

                LogLike(n)=0;
                for w=1:units
%                     SPET=Classifier(1).TotalRastersFm2(w).Units(m).Trials(ival(m,w)).spet;
                    SPET=Classifier(1).TotalRastersFm2fifty(w).Units(m).Trials(ival(m,w)).spet;
                    SPET=SPET(SPET~=0);
                    LogLike(n)=LogLike(n)+sum(log10(PSTHi(n,ceil(SPET/Fs*1000),w)));

                end
            end
                [hit,index]=max(real(LogLike));
                Results(m,index,l)=1;
end
end

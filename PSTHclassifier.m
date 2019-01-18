% function [Results]=PSTHclassifier(Rasterall,FsHist,NB)
%   
%       DESCRIPTION     : This function uses a bernoulli bayesian calssifier to estimate the percent correct of single 
%                         trails compared to PSTH templates. 
%       Rasterall       : All rasters from one neuron. Rastergram Data Structure
%                         spet: spike event time 
%                         Fs: sampling rate
%       FsHist          : Desired Sampling Rate for cycle PSTH (Hz)
%       NB              : Number of bootstraps to be used. 
%
%Returned Values
%       Results         : Results from comparing single trials to the PSTH
%                         templates. (input sound,outputsound,bootstrap
%                         number) 
% (C) Ahmad Osman & Monty Escabi, Jan 2017

function [Results]=PSTHclassifier(Rasterall,FsHist,NB)
%Initialize variables
N=length(Rasterall); %Number of sounds
L=length(Rasterall(1).RR); %Number of trials 
Results=zeros(N,N,NB);  %Results matrix containig 0 (miss) or 1 (hit)
Fs=Rasterall(1).RR(1).Fs;   %spike train sampling rate 


for l=1:NB %Bootstraps
        %Randomly Find model and validation trials for all sounds (for
        %Bootstrap). Also generates models for validation
        for n=1:N   %output sounds
            %Find model and validation trials
            imod(n,:)=randsample(L,L-1)';   %randomly remove 1 trial - model trials
            ival(n)=sum(1:L)-sum(imod(n,:)); %validation trial
            
            %Generate Model Cycle Raster
            [RASTER]=raster2cycleraster(Rasterall(n).RR(imod(n,:)),2,1,0,0);% creating cycle raster 2 Hz Fm 
            [PSTHi(n,:)]=cyclepsthmodel(RASTER,FsHist,1000,16); %creating PSTH template
            PSTHi(n,:)=PSTHi(n,:)/sum(PSTHi(n,:)); %normalizing 
        end
        
        %Validating for current bootstrap itteration
        for m=1:N   %input sounds
            LogLike=zeros(1,N);
            for n=1:N   %ouput sounds
                LogLike(n)=sum(log10(PSTHi(n,ceil(Rasterall(m).RR(ival(m)).spet/Fs*1000))));
            end
            [hit,index]=max(real(LogLike));
            Results(m,index,l)=1;
        end
end




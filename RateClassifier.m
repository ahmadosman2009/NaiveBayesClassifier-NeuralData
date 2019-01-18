% function [Results] = RateClassifier(Rasterall,NB,Fm)
%   
%       DESCRIPTION     : This function uses a bernoulli bayesian calssifier to estimate the percent correct of single 
%                         trails compared to cPSTH templates using rate code. 
%       Rasterall       : All rasters from one neuron. Rastergram Data Structure
%                         spet: spike event time 
%                         Fs: sampling rate
%       NB              : Number of bootstraps to be used. 
%       Fm              : Modulation frequency (Hz)
%
%Returned Values
%       Results         : Results from comparing single trials to the PSTH
%                         templates. (input sound,outputsound,bootstrap
%                         number) 
% (C) Ahmad Osman & Monty Escabi, Feb 2017

function [Results] = RateClassifier(Rasterall,NB,Fm)

N=length(Rasterall); %Number of sounds
L=length(Rasterall(1).RR); %Number of trials 
Results=zeros(N,N,NB);  %Results matrix containig 0 (miss) or 1 (hit)
Fs=Rasterall(1).RR(1).Fs;   %spike train sampling rate
Tp=1/Fm;    %Sound Period

for l=1:NB %Bootstraps
    %Generating Model
    for n=1:N   %sounds
                %Find model and validation trials
                imod(n,:)=randsample(L,L-1)';   %randomly remove 1 trial - model trials
                ival(n)=sum(1:L)-sum(imod(n,:)); %validation trial

                %Genearging raster for model trials only
                [RAS]=rasterexpand(Rasterall(n).RR(imod(n,:)),Fm,4);
                lambda(n)= mean(mean(RAS(:,1:4)));                          %Average firing rate for model trials - this is our model (Poisson)
    end

    %Classifying
    for m=1:N   %input sounds
%             LogLike=zeros(1,N);
            for n=1:N   %ouput sounds (model)
                Loglike(n)=0;
                 %Validating for current bootstrap itteration
                [RASv]=rasterexpand(Rasterall(m).RR(ival(m)),Fm,4);     %Validatoin response - single trial with multiple cycles
                SpikeCountV=RASv(:,1:4)*Tp;                             %Spike count for each validation cycle
                P=poisspdf(SpikeCountV,lambda(n)*Tp);
                LogLike(n)= sum(log10(P)); 
            end
            
    %Makind decision to maximize log-like
    [hit,index]=max(real(LogLike));
    Results(m,index,l)=1;
    end  
end 

%function [PSTHi,newFsd,taxis]=cyclepsthmodel(RASTER,FsHist,Fsd,L)
%           [PSTHi,newFsd,taxis]=cyclepsthmodel(RASTER,200,1000,4)
%       DESCRIPTION     : Generates an interpolated PSTH from a RASTER
%       RASTER          : Rastergram Data Structure
%                         spet: spike event time 
%                         Fs: sampling rate
%       FsHist          : Desired Sampling Rate for cycle PSTH (Hz)
%       Fsd             : Sampling rate for interplolation (Hz). 
%       L               : Number of desired interplolated PSTH cycles 
%Returned Values
%
%       PSTHi           : Interpolated PSTH
%       newFsd          : new sampling rate rounded up (ceil) to create
%                         integer number of samples 
%       taxis           : New time axis for final PSTH (sec). 
% (C) Ahmad Osman & Monty Escabi, Dec 2016
%
function [PSTHi,newFsd,taxis]=cyclepsthmodel(RASTER,FsHist,Fsd,L)
T=RASTER(1).T; %duration of cycle 
[RAS]=rasterexpand(RASTER,FsHist,T);
PSTHc=mean(RAS)/FsHist; %1 cycle PSTH
PSTHc=[PSTHc PSTHc PSTHc]; %3 cycle PSTH
d=(0:length(PSTHc)-1)/FsHist;   %histogram time axis

newFsd=ceil(Fsd*T)/T;  %Check and change sampling rate so there are an integer number of samples in one cycle
N=round(newFsd*T*3);   %Number of samples for 3 cycles
t=(0:N-1)/newFsd;    %Intermpolation time axis
PSTHcInterp = interp1(d,PSTHc,t,'spline');
R=length(PSTHcInterp)/3; %lentgh of one cycle 
PSTHi=PSTHcInterp(R:2*R-1); %selecting the middle cycle 
PSTHi=repmat(PSTHi,1,L); 
taxis=0:1/newFsd:length(PSTHi)/newFsd-1/newFsd; 

%%Convoloving with gaussian to remove zero terms 
NN=length(PSTHi)/2;
n=-NN:NN-1;
G=exp(-n.^2/2/(2)^2); 
G=G/sum(G);
R=xcorrcircular(PSTHi,G);
PSTHi=fftshift(R); 




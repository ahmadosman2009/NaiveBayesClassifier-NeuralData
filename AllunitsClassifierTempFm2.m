%%%%%%%%%%%%%%%%%%%%%%%%% Extracting units from acstats then Running Temp Classifier  %%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%                   Loading acstats                           %%%%%%%%%%%%%%%%%%%%%%%%% 
% load('C:\Users\The Shit\Desktop\Lab\Discrimination\acstatsv2ss_F2.mat')
load('E:\PhD Lab\Discrimination\acstatsv2ss_F2.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%                Loading in Variable                          %%%%%%%%%%%%%%%%%%%%%%%%% 
[FM, FC, FMAxis, FCAxis, Nbw, Nfm] = makeFMAxis('splinemodnoise64Hz_param.mat');    

for i = 1:numel(acstats)
    acstats(i).number = i;
end

in(1).stats = [];
region = cell2mat(s2a(acstats, 'registry(1,12)'));
carrier = cell2mat(s2a(acstats, 'registry(1,7)'));
CF = cell2mat(s2a(acstats, 'registry(1,15)'));
BF = cell2mat(s2a(acstats, 'registry(1,14)'));

in(2).stats = acstats(1,region == 1 & carrier == 2 & CF > 4000 & CF < 250000 & abs(s2a(acstats, 'spikeWidth2')) > 0.2 & abs(s2a(acstats, 'spikeWidth2')) < 0.7); %A1 noise
in(3).stats = acstats(1,region == 2 & carrier == 2 & CF > 4000 & CF < 25000 & abs(s2a(acstats, 'spikeWidth2')) > 0.2 & abs(s2a(acstats, 'spikeWidth2')) < 0.7); %VAF noise
in(4).stats = acstats(1,region == 3 & carrier == 2 & CF > 4000 & CF < 250000 & abs(s2a(acstats, 'spikeWidth2')) > 0.2 & abs(s2a(acstats, 'spikeWidth2')) < 0.7); %SRAF noise

%%%%%%%%%%%%%%%%%%%%%%%%%          Saving 2Hz Fm Conditions into Classifier Var       %%%%%%%%%%%%%%%%%%%%%%%%% 
m=length(in(2).stats); %A1 units 
for n=1:m
    for w=1:10
        Classifier(1).RastersFm2(n).Units(w).Trials=in(2).stats(n).stim(1,w).RASTER;
    end
end 

m=length(in(3).stats); %VAF untis 
for n=1:m
    for w=1:10
        Classifier(2).RastersFm2(n).Units(w).Trials=in(3).stats(n).stim(1,w).RASTER;
    end
end 

m=length(in(4).stats); %SRAF units 
for n=1:m
    for w=1:10
        Classifier(3).RastersFm2(n).Units(w).Trials=in(4).stats(n).stim(1,w).RASTER;
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%          Running SU Classifier              %%%%%%%%%%%%%%%%%%%%%%%%% 

%%% calc of every unit in each region for temporal code 
FsHist=200;
NB=200;
L=length(Classifier(1).RastersFm2); %A1 units 
N=length(Classifier(1).RastersFm2(1).Units); %10 sounds 
for k=1:L
    
    for p=1:N %10 sounds
    Rasterall(p).RR=Classifier(1).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(1).TempSU(k).Unit]=PSTHclassifier(Rasterall,FsHist,NB);
    k
end
 
%VAF units 
L=length(Classifier(2).RastersFm2);  %VAF units 
N=length(Classifier(2).RastersFm2(1).Units); %10 sounds 
for k=1:L
    
    for p=1:N %10 sounds
    Rasterall(p).RR=Classifier(2).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(2).TempSU(k).Unit]=PSTHclassifier(Rasterall,FsHist,NB);
    k
end
    
%SRAF units 
L=length(Classifier(3).RastersFm2); %SRAF units 
N=length(Classifier(3).RastersFm2(1).Units); %10 sounds 
for k=1:L
    
    for p=1:N %10 sounds
    Rasterall(p).RR=Classifier(3).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(3).TempSU(k).Unit]=PSTHclassifier(Rasterall,FsHist,NB);
    k
end


%%%%%%%%%%%%%%%%%%%%%%%%%                   Mean of the bootstraps           %%%%%%%%%%%%%%%%%%%%%%%%% 
L=length(Classifier(1).RastersFm2);
for g=1:L
Classifier(1).TempSU(g).sum=mean(Classifier(1).TempSU(g).Unit,3);
end 

L=length(Classifier(2).RastersFm2);
for g=1:L
Classifier(2).TempSU(g).sum=mean(Classifier(2).TempSU(g).Unit,3);
end 


L=length(Classifier(3).RastersFm2);
for g=1:L
Classifier(3).TempSU(g).sum=mean(Classifier(3).TempSU(g).Unit,3);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%               SU avg per region                    %%%%%%%%%%%%%%%%%%%%%%%%% 
L=length(Classifier(1).RastersFm2);
Classifier(1).TempSU(1).total=0;
for s=1:L
    Classifier(1).TempSU(1).total=Classifier(1).TempSU(1).total+Classifier(1).TempSU(s).sum;
end
    Classifier(1).TempSU(1).total=Classifier(1).TempSU(1).total/L;
 
    
    
%   
%     
%     
%     
%     
%     
%     
%    
%     
%     
%     

    

    L=length(Classifier(2).RastersFm2);
Classifier(2).TempSU(1).total=0;
for s=1:L
    Classifier(2).TempSU(1).total=Classifier(2).TempSU(1).total+Classifier(2).TempSU(s).sum;
end
    Classifier(2).TempSU(1).total=Classifier(2).TempSU(1).total/L;
    

      L=length(Classifier(3).RastersFm2);
Classifier(3).TempSU(1).total=0;
for s=1:L
    Classifier(3).TempSU(1).total=Classifier(3).TempSU(1).total+Classifier(3).TempSU(s).sum;
end
    Classifier(3).TempSU(1).total=Classifier(3).TempSU(1).total/L;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%            Running Pop avg per region Classifier                       %%%%%%%%%%%%%%%%%%%%%%%%%    


for n=1:50
[Classifier(1).TempPop(n).Matrix]=TemporalPopulationClassifier(Classifier,1,200,n);
Classifier(1).TempPop(n).Mean=mean(Classifier(1).TempPop(n).Matrix,3);
Classifier(1).TempPop(n).DiagMean=mean(diag(Classifier(1).TempPop(n).Mean));
n 
%total 2 hrs
end

for n=1:50
[Classifier(2).TempPop(n).Matrix]=TemporalPopulationClassifierFm2(Classifier,2,200,n);
Classifier(2).TempPop(n).Mean=mean(Classifier(2).TempPop(n).Matrix,3);
Classifier(2).TempPop(n).DiagMean=mean(diag(Classifier(2).TempPop(n).Mean));
n
end

for n=1:50
[Classifier(3).TempPop(n).Matrix]=TemporalPopulationClassifierFm2(Classifier,3,200,n);
Classifier(3).TempPop(n).Mean=mean(Classifier(3).TempPop(n).Matrix,3);
Classifier(3).TempPop(n).DiagMean=mean(diag(Classifier(3).TempPop(n).Mean));
n
end

L=50;
    Classifier(1).TempPop(1).total=0;
for s=1:L
    Classifier(1).TempPop(1).total=Classifier(1).TempPop(1).total+Classifier(1).TempPop(s).Mean;
end
    Classifier(1).TempPop(1).total=Classifier(1).TempPop(1).total/L;

    Classifier(2).TempPop(1).total=0;
for s=1:L
    Classifier(2).TempPop(1).total=Classifier(2).TempPop(1).total+Classifier(2).TempPop(s).Mean;
end
    Classifier(2).TempPop(1).total=Classifier(2).TempPop(1).total/L;
    
    Classifier(3).TempPop(1).total=0;
for s=1:L
    Classifier(3).TempPop(1).total=Classifier(3).TempPop(1).total+Classifier(3).TempPop(s).Mean;
end
    Classifier(3).TempPop(1).total=Classifier(3).TempPop(1).total/L;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%             Arranging Data For Total Pop              %%%%%%%%%%%%%%%%%%%%%%%%%        
    
for g=1:50
Classifier(1).TotalRastersFm2(g).Units=Classifier(1).RastersFm2(g).Units
end 

for g=1:50
Classifier(1).TotalRastersFm2(g+50).Units=Classifier(2).RastersFm2(g).Units
end

for g=1:50
Classifier(1).TotalRastersFm2(g+100).Units=Classifier(3).RastersFm2(g).Units
end     
    

%%%%%%%%%%%%%%%%%%%%%%%%%            Running Total Pop avg Classifier            %%%%%%%%%%%%%%%%%%%%%%%%%   
for n=1:150
[Classifier(1).TotalPopulation(n).Matrix]=TemporalTotalPopulationClassifierFm2(Classifier,200,n);
Classifier(1).TotalPopulation(n).Mean=mean(Classifier(1).TotalPopulation(n).Matrix,3);
Classifier(1).TotalPopulation(n).DiagMean=mean(diag(Classifier(1).TotalPopulation(n).Mean));
n
%start 10:37 am
end


%%%% re-running the total pop for n=50 instead of n=150 
uu=randsample(50,20);
zz=randsample(50,15)+50;
qq=randsample(50,15)+100;
randomfifty=[uu' zz' qq'];

for v=1:50
Classifier(1).TotalRastersFm2fifty(v).Units=Classifier(1).TotalRastersFm2(randomfifty(v)).Units;
end 

for n=1:50
[Classifier(1).TotalPopulationfifty(n).Matrix]=TemporalTotalPopulationClassifierFm2(Classifier,200,n);
Classifier(1).TotalPopulationfifty(n).Mean=mean(Classifier(1).TotalPopulationfifty(n).Matrix,3);
Classifier(1).TotalPopulationfifty(n).DiagMean=mean(diag(Classifier(1).TotalPopulationfifty(n).Mean));
n
%start 5:40 pm
end

L=50;
    Classifier(1).TotalPopulationfifty(1).total=0;
for s=1:L
    Classifier(1).TotalPopulationfifty(1).total=Classifier(1).TotalPopulationfifty(1).total+Classifier(1).TotalPopulationfifty(s).Mean;
end
    Classifier(1).TotalPopulationfifty(1).total=Classifier(1).TotalPopulationfifty(1).total/L;
    
imagesc(Classifier(1).TotalPopulationfifty(1).total,[0 1]),colorbar,set(gca,'YDir','normal')



L=150;
    Classifier(1).TotalPopulation(1).total=0;
for s=1:L
    Classifier(1).TotalPopulation(1).total=Classifier(1).TotalPopulation(1).total+Classifier(1).TotalPopulation(s).Mean;
end
    Classifier(1).TotalPopulation(1).total=Classifier(1).TotalPopulation(1).total/L;
    
imagesc(Classifier(1).TotalPopulation(1).total,[0 1]),colorbar,set(gca,'YDir','normal')

 %%%% Figure 8-I: population temporal code A1
imagesc(Classifier(1).TempPop(1).total,[0 1]),colorbar,set(gca,'YDir','normal')
title('A1 50 units')    

 %%%% Figure 8-J: population temporal code VAF
imagesc(Classifier(2).TempPop(1).total,[0 1]),colorbar,set(gca,'YDir','normal')
title('VAF 50 units')

 %%%% Figure 8-K: population temporal code SRAF
imagesc(Classifier(3).TempPop(1).total,[0 1]),colorbar,set(gca,'YDir','normal')
title('SRAF 50 units')

 %%%% Figure 8-A: single neuron temporal code A1
figure       
imagesc(Classifier(1).TempSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')  

 %%%% Figure 8-B: single neuron temporal code VAF
figure
imagesc(Classifier(2).TempSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal') 

 %%%% Figure 8-C: single neuron temporal code SRAF
figure
imagesc(Classifier(3).TempSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')  
    
for h=1:54
imagesc(Classifier(1).TempSU(h).sum2,[0 0.5]),colorbar,set(gca,'YDir','normal')
title(['cell number' num2str(h)])
pause
end

for t=1:54
Classifier(1).TempSU(t).Unit=Classifier(1).TempSU(t).Unit2;
end
    

for k=1:54
A1diagmeanSUsFm2(k,:)=diag(Classifier(1).TempSU(k).sum);
end
A1stdSUsFm2=bootstrp(1000,'mean',A1diagmeanSUsFm2);
figure
errorbar(1:10,mean(A1diagmeanSUsFm2),std(A1stdSUsFm2))
xlabel('Sound shape')
ylabel('% correct')
title('A1 population % correct')
ylim([0 1])

for k=1:104
VAFdiagmeanSUsFm2(k,:)=diag(Classifier(2).TempSU(k).sum);
end
VAFstdSUsFm2=bootstrp(1000,'mean',VAFdiagmeanSUsFm2);
errorbar(1:10,mean(VAFdiagmeanSUsFm2),std(VAFstdSUsFm2))
xlabel('Sound shape')
ylabel('% correct')
title('VAF population % correct')
ylim([0 1])

for k=1:65
SRAFdiagmeanSUsFm2(k,:)=diag(Classifier(3).TempSU(k).sum);
end
SRAFstdSUsFm2=bootstrp(1000,'mean',SRAFdiagmeanSUsFm2);
errorbar(1:10,mean(SRAFdiagmeanSUsFm2),std(SRAFstdSUsFm2))
xlabel('Sound shape')
ylabel('% correct')
title('SRAF population % correct')
ylim([0 1])

 %%%% Figure 8-D: single neuron temporal code A1,VAF, and SRAF matrix
 %%%% diagonal
figure
Fc=[2,4,6,8,11,16,23,32,45,64];
errorbar(Fc,(mean(A1diagmeanSUsFm2)).*100,(std(A1stdSUsFm2)).*100,'b','LineWidth', 2)
hold on 
errorbar(Fc,(mean(VAFdiagmeanSUsFm2)).*100,(std(VAFstdSUsFm2)).*100,'g','LineWidth', 2)
errorbar(Fc,(mean(SRAFdiagmeanSUsFm2)).*100,(std(SRAFstdSUsFm2)).*100,'r','LineWidth', 2)
xlabel('Fc (Hz)')
ylabel('% correct')
title('A1, VAF, SRAF total SU % correct Temporal Code')
ylim([0 100])
legend('A1', 'VAF', 'SRAF')
set(gca,'xscale', 'log')

imagesc(Classifier(1).RateSU(1).total,[0 0.6]),colorbar,set(gca,'YDir','normal')

figure
imagesc(Classifier(2).RateSU(1).total,[0 0.6]),colorbar,set(gca,'YDir','normal')
 
figure
imagesc(Classifier(3).RateSU(1).total,[0 0.6]),colorbar,set(gca,'YDir','normal')

Fc=[2,4,6,8,11,16,23,32,45,64];
for k=1:50
A1diagmeanFm2(k,:)=diag(Classifier(1).TempPop(k).Mean);
end
A1diagbootFm2=bootstrp(1000,'mean',A1diagmeanFm2);
errorbar(Fc,(mean(A1diagmeanFm2)).*100,(std(A1diagbootFm2)).*100)
xlabel('Fc (Hz)')
ylabel('% correct')
title('A1 population % correct')
ylim([0 100])


for k=1:50
VAFdiagmeanFm2(k,:)=diag(Classifier(2).TempPop(k).Mean);
end
VAFdiagbootFm2=bootstrp(1000,'mean',VAFdiagmeanFm2);
errorbar(Fc,(mean(VAFdiagmeanFm2)).*100,(std(VAFdiagbootFm2)).*100)
xlabel('Fc (Hz)')
ylabel('% correct')
title('VAF population % correct')
ylim([0 100])


for k=1:50
SRAFdiagmeanFm2(k,:)=diag(Classifier(3).TempPop(k).Mean);
end
SRAFdiagbootFm2=bootstrp(1000,'mean',SRAFdiagmeanFm2);
errorbar(Fc,(mean(SRAFdiagmeanFm2)).*100,(std(SRAFdiagbootFm2)).*100)
% xlabel('Sound shape')
xlabel('Fc (Hz)')
ylabel('% correct')
title('SRAF population % correct')
ylim([0 100])

%%%% Figure 8-L (part1): population temporal code A1,VAF,SRAF diagnoal
errorbar(Fc,(mean(A1diagmeanFm2)).*100,(std(A1diagbootFm2)).*100,'b','LineWidth', 2)
hold on 
errorbar(Fc,(mean(VAFdiagmeanFm2)).*100,(std(VAFdiagbootFm2)).*100,'g','LineWidth', 2)
errorbar(Fc,(mean(SRAFdiagmeanFm2)).*100,(std(SRAFdiagbootFm2)).*100,'r','LineWidth', 2)
xlabel('Fc (Hz)')
ylabel('% correct')
title('A1, VAF, SRAF population % correct Temporal Code')
ylim([0 100])
legend('A1', 'VAF', 'SRAF')
set(gca,'xscale', 'log')



%%%% Figure 8-L (part1+part2): population temporal+rate code A1,VAF,SRAF diagnoal
errorbar(Fc,(mean(A1diagmeanFm2)).*100,(std(A1diagbootFm2)).*100,'b','LineWidth', 2)
hold on 
errorbar(Fc,(mean(VAFdiagmeanFm2)).*100,(std(VAFdiagbootFm2)).*100,'g','LineWidth', 2)
errorbar(Fc,(mean(SRAFdiagmeanFm2)).*100,(std(SRAFdiagbootFm2)).*100,'r','LineWidth', 2)
errorbar(Fc,mean(A1diagmeanRFm2).*100,std(A1diagbootR).*100,'-ob','LineWidth', 2)
errorbar(Fc,mean(VAFdiagmeanRFm2).*100,std(VAFdiagbootRFm2).*100,'-og','LineWidth', 2)
errorbar(Fc,mean(SRAFdiagmeanRFm2).*100,std(SRAFdiagbootRFm2).*100,'-or','LineWidth', 2)
xlabel('Fc (Hz)')
ylabel('% correct')
title('A1, VAF, SRAF population % correct Temporal & Rate Code')
ylim([0 100])
legend('A1 Temporal', 'VAF Temporal', 'SRAF Temporal','A1 Rate', 'VAF Rate', 'SRAF Rate')
set(gca,'xscale', 'log')


%%% plotting 50 A1 neurons vs %correct with error shade 
for k=1:50
A1diagFm2(k,:)=diag(Classifier(1).TempPop(k).Mean);
A1bootFm2(k,:)=bootstrp(1000,'mean',A1diagFm2(k,:));
A1eFm2(k)=std(A1bootFm2(k,:));
end
x=1:50;
x=x';
A1yFm2=mean(A1diagFm2,2)*100;
A1eFm2=A1eFm2'*100;
shadedErrorBar(x,A1yFm2,A1eFm2,'b',.5)

[A1ysortFm2,IsortFm2Fm2]=sort(A1yFm2,1);
for h=1:50
A1esortFm2(h)=A1eFm2(IsortFm2Fm2(h));
end
shadedErrorBar(x,A1ysortFm2,A1esortFm2,'b')

title('A1')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])

%%% plotting 50 VAF neurons vs %correct with error shade 
for k=1:50
VAFdiagFm2(k,:)=diag(Classifier(2).TempPop(k).Mean);
VAFbootFm2(k,:)=bootstrp(1000,'mean',VAFdiagFm2(k,:));
VAFeFm2(k)=std(VAFbootFm2(k,:));
end
x=1:50;
x=x';
VAFyFm2=mean(VAFdiagFm2,2)*100;
VAFeFm2=VAFeFm2'*100;
shadedErrorBar(x,VAFyFm2,VAFeFm2,'g',.5)

[VAFysortFm2,IsortFm2Fm2]=sort(VAFyFm2,1);
for h=1:50
VAFesortFm2(h)=VAFeFm2(IsortFm2Fm2(h));
end
shadedErrorBar(x,VAFysortFm2,VAFesortFm2,'g')

title('VAF')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])


%%% plotting 50 SRAF neurons vs %correct with error shade 
for k=1:50
SRAFdiagFm2(k,:)=diag(Classifier(3).TempPop(k).Mean);
SRAFbootFm2(k,:)=bootstrp(1000,'mean',SRAFdiagFm2(k,:));
SRAFeFm2(k)=std(SRAFbootFm2(k,:));
end
x=1:50;
x=x';
SRAFyFm2=mean(SRAFdiagFm2,2)*100;
SRAFeFm2=SRAFeFm2'*100;
shadedErrorBar(x,SRAFyFm2,SRAFeFm2,'r',.5)

[SRAFysortFm2,IsortFm2Fm2]=sort(SRAFyFm2,1);
for h=1:50
SRAFesortFm2(h)=SRAFeFm2(IsortFm2Fm2(h));
end
shadedErrorBar(x,SRAFysortFm2,SRAFesortFm2,'r',.5)

title('SRAF')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])

%%%% Figure 8-P (part1): population temporal # of neurons
shadedErrorBar(x,A1ysortFm2,A1esortFm2,'b')
hold on 
shadedErrorBar(x,VAFysortFm2,VAFesortFm2,'g')
shadedErrorBar(x,SRAFysortFm2,SRAFesortFm2,'r')
title('A1,VAF,SRAF')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])

%%%% Figure 8-P: population temporal and rate code # of
%%%% neurons vs classification 
shadedErrorBar(x,A1ysortFm2,A1esortFm2,'b')
hold on 
shadedErrorBar(x,VAFysortFm2,VAFesortFm2,'g')
shadedErrorBar(x,SRAFysortFm2,SRAFesortFm2,'r')
shadedErrorBar(x,A1ysortRFm2,A1esortRFm2,':b')
shadedErrorBar(x,VAFysortRFm2,VAFesortRFm2,':g')
shadedErrorBar(x,SRAFysortRFm2,SRAFesortRFm2,':r')
title('A1,VAF,SRAF Temporal and Rate code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])




Fc=[2,4,6,8,11,16,23,32,45,64];
for k=1:150
AlldiagmeanTFm2(k,:)=diag(Classifier(1).TotalPopulation(k).Mean);
end
AlldiagbootTFm2=bootstrp(1000,'mean',AlldiagmeanTFm2);
errorbar(Fc,(mean(AlldiagmeanTFm2)).*100,(std(AlldiagbootTFm2)).*100,'r','LineWidth', 2)
% xlabel('Sound shape')
xlabel('Sound Shape Fc (Hz)')
ylabel('% correct')
title('All population % correct Temporal')
ylim([0 100])
set(gca,'xscale', 'log')

errorbar(Fc,(mean(AlldiagmeanTFm2)).*100,(std(AlldiagbootTFm2)).*100,'r','LineWidth', 2)
hold on 
errorbar(Fc,(mean(AlldiagmeanRFm2)).*100,(std(AlldiagbootRFm2)).*100,'b*-','LineWidth', 2)
xlabel('Sound Shape Fc (Hz)')
ylabel('% correct')
title('All population % correct Rate and Temporal')
ylim([0 100])
set(gca,'xscale', 'log')



for k=1:150
AlldiagFm2(k,:)=diag(Classifier(1).TotalPopulation(k).Mean);
AllbootFm2(k,:)=bootstrp(1000,'mean',AlldiagFm2(k,:));
AlleFm2(k)=std(AllbootFm2(k,:));
end
x=1:150;
x=x';
AllyFm2=mean(AlldiagFm2,2)*100;
AlleFm2=AlleFm2'*100;
shadedErrorBar(x,AllyFm2,AlleFm2,'r',.5)

[AllysortFm2,IsortFm2]=sort(AllyFm2,1);
for h=1:150
AllesortFm2(h)=AlleFm2(IsortFm2(h));
end
shadedErrorBar(x,AllysortFm2,AllesortFm2,'r',.5)

title('A1,VAF,SRAF Temporal Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])



shadedErrorBar(x,AllysortFm2,AllesortFm2,'r')
hold on 
shadedErrorBar(x,AllysortRFm2,AllesortRFm2,'*-b')
title('A1,VAF,SRAF Temporal and Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])

    
NB=200;
L=length(Classifier(1).RastersFm2); %A1 units 
N=length(Classifier(1).RastersFm2(1).Units); %10 sounds 
Fm=2;
for k=1:L
    
    for p=1:N %9 sounds
    Rasterall(p).RR=Classifier(1).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(1).RateSU(k).Unit]= RateClassifier(Rasterall,NB,Fm);
    k
end
 
%%%%%%%%%%% 
%VAF units 
L=length(Classifier(2).RastersFm2); %A1 units 
N=length(Classifier(2).RastersFm2(1).Units); %10 sounds 
for k=1:L
    
    for p=1:N %9 sounds
    Rasterall(p).RR=Classifier(2).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(2).RateSU(k).Unit]= RateClassifier(Rasterall,NB,Fm)
    k
end
    
%SRAF units 
L=length(Classifier(3).RastersFm2); %A1 units 
N=length(Classifier(3).RastersFm2(1).Units); %10 sounds 
for k=1:L
    
    for p=1:N %9 sounds
    Rasterall(p).RR=Classifier(3).RastersFm2(k).Units(p).Trials;
    end
    [Classifier(3).RateSU(k).Unit]= RateClassifier(Rasterall,NB,Fm)
    k
end


%mean of the 200 bootstraps 
L=length(Classifier(1).RastersFm2);
for g=1:L
Classifier(1).RateSU(g).sum=mean(Classifier(1).RateSU(g).Unit,3);
end 

L=length(Classifier(2).RastersFm2);
for g=1:L
Classifier(2).RateSU(g).sum=mean(Classifier(2).RateSU(g).Unit,3);
end 

L=length(Classifier(3).RastersFm2);
for g=1:L
Classifier(3).RateSU(g).sum=mean(Classifier(3).RateSU(g).Unit,3);
end 

 %%%% mean of the whole region 
L=length(Classifier(1).RastersFm2);
Classifier(1).RateSU(1).total=0;
for s=1:L
    Classifier(1).RateSU(1).total=Classifier(1).RateSU(1).total+Classifier(1).RateSU(s).sum;
end
    Classifier(1).RateSU(1).total=Classifier(1).RateSU(1).total/L;
    
    L=length(Classifier(2).RastersFm2);
Classifier(2).RateSU(1).total=0;
for s=1:L
    Classifier(2).RateSU(1).total=Classifier(2).RateSU(1).total+Classifier(2).RateSU(s).sum;
end
    Classifier(2).RateSU(1).total=Classifier(2).RateSU(1).total/L;
    
    L=length(Classifier(3).RastersFm2);
Classifier(3).RateSU(1).total=0;
for s=1:L
    Classifier(3).RateSU(1).total=Classifier(3).RateSU(1).total+Classifier(3).RateSU(s).sum;
end
    Classifier(3).RateSU(1).total=Classifier(3).RateSU(1).total/L;
    
     %%%% Figure 8-E: single neuron rate code A1
    figure
    imagesc(Classifier(1).RateSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')
    
     %%%% Figure 8-F: single neuron rate code VAF
    figure
    imagesc(Classifier(2).RateSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')
    
    %%%% Figure 8-G: single neuron rate code SRAF
    figure
    imagesc(Classifier(3).RateSU(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')
    
    
    
    
    
    
    
    
    
    
     
   for k=1:54
A1diagmeanSUsRFm2(k,:)=diag(Classifier(1).RateSU(k).sum);
   end
A1stdSUsRFm2=bootstrp(1000,'mean',A1diagmeanSUsRFm2);
figure
errorbar(1:10,mean(A1diagmeanSUsRFm2),std(A1stdSUsRFm2))
xlabel('Sound shape')
  ylabel('% correct')
  title('A1 population % correct Rate Code')
  ylim([0 0.5])
   
   
   
   for k=1:104
VAFdiagmeanSUsRFm2(k,:)=diag(Classifier(2).RateSU(k).sum);
   end
VAFstdSUsRFm2=bootstrp(1000,'mean',VAFdiagmeanSUsRFm2);
figure
errorbar(1:10,mean(VAFdiagmeanSUsRFm2),std(VAFstdSUsRFm2))
xlabel('Sound shape')
  ylabel('% correct')
  title('VAF population % correct Rate Code')
  ylim([0 0.5])
   
   for k=1:65
SRAFdiagmeanSUsRFm2(k,:)=diag(Classifier(3).RateSU(k).sum);
   end
SRAFstdSUsRFm2=bootstrp(1000,'mean',SRAFdiagmeanSUsRFm2);
figure
errorbar(1:10,mean(SRAFdiagmeanSUsRFm2),std(SRAFstdSUsRFm2))
xlabel('Sound shape')
  ylabel('% correct')
  title('SRAF population % correct Rate Code')
  ylim([0 0.5])
   
   %%%% Figure 8-H: single neuron rate code A1, VAF, SRAF diagonal
 Fc=[2,4,6,8,11,16,23,32,45,64];
errorbar(Fc,mean(A1diagmeanSUsRFm2)*100,std(A1stdSUsRFm2)*100, '*-', 'LineWidth', 2)
hold on 
errorbar(Fc,mean(VAFdiagmeanSUsRFm2)*100,std(VAFstdSUsRFm2)*100, 'g*-', 'LineWidth', 2)
errorbar(Fc,mean(SRAFdiagmeanSUsRFm2)*100,std(SRAFstdSUsRFm2)*100, 'r*-', 'LineWidth', 2)
xlabel('Fc (Hz)')
  ylabel('% correct')
  title('A1,VAF,SRAF population % correct Rate Code')
  ylim([0 50])
   legend('A1','VAF','SRAF')
   set(gca,'xscale', 'log')

  
    
        
for n=1:50
[Classifier(1).RatePop(n).Matrix]=RatePopulationClassifier3Fm2Hz(Classifier,1,200,2,n);
Classifier(1).RatePop(n).Mean=mean(Classifier(1).RatePop(n).Matrix,3); 
Classifier(1).RatePop(n).DiagMean=mean(diag(Classifier(1).RatePop(n).Mean));
n 
%3:47 pm
%4:14 pm 25 done
%5:43 pm 50 
%total 2 hrs
end
    
 L=50;
Classifier(1).RatePop(1).total=0;
for s=1:L
    Classifier(1).RatePop(1).total=Classifier(1).RatePop(1).total+Classifier(1).RatePop(s).Mean;
end
    Classifier(1).RatePop(1).total=Classifier(1).RatePop(1).total/L;
    
    for n=1:50
    tic
[Classifier(2).RatePop(n).Matrix2]=RatePopulationClassifier3Fm2Hz(Classifier,2,200,2,n);
Classifier(2).RatePop(n).Mean2=mean(Classifier(2).RatePop(n).Matrix2,3);
Classifier(2).RatePop(n).DiagMean2=mean(diag(Classifier(2).RatePop(n).Mean2));
toc 
n 
    end   
    
for n=1:50
    tic
[Classifier(3).RatePop(n).Matrix]=RatePopulationClassifier3Fm2Hz(Classifier,3,200,2,n);
Classifier(3).RatePop(n).Mean=mean(Classifier(3).RatePop(n).Matrix,3);
Classifier(3).RatePop(n).DiagMean=mean(diag(Classifier(3).RatePop(n).Mean));
toc 
n 
end


for n=1:150
tic
[Classifier(1).TotalPopulationRate(n).Matrix]=RateTotalPopulationClassifier3Fm2Hz(Classifier,200,2,n);
Classifier(1).TotalPopulationRate(n).Mean=mean(Classifier(1).TotalPopulationRate(n).Matrix,3);
Classifier(1).TotalPopulationRate(n).DiagMean=mean(diag(Classifier(1).TotalPopulationRate(n).Mean));
toc 
n
end

%%%% re-running the total pop for n=50 instead of n=150 
for n=1:50
[Classifier(1).TotalPopulationRatefifty(n).Matrix]=RateTotalPopulationClassifier3Fm2Hz(Classifier,200,2,n);
Classifier(1).TotalPopulationRatefifty(n).Mean=mean(Classifier(1).TotalPopulationRatefifty(n).Matrix,3);
Classifier(1).TotalPopulationRatefifty(n).DiagMean=mean(diag(Classifier(1).TotalPopulationRatefifty(n).Mean));
n
%start 11:07 am
%end  12:50 pm 
end

L=50;
    Classifier(1).TotalPopulationRatefifty(1).total=0;
for s=1:L
    Classifier(1).TotalPopulationRatefifty(1).total=Classifier(1).TotalPopulationRatefifty(1).total+Classifier(1).TotalPopulationRatefifty(s).Mean;
end
    Classifier(1).TotalPopulationRatefifty(1).total=Classifier(1).TotalPopulationRatefifty(1).total/L;
    
imagesc(Classifier(1).TotalPopulationRatefifty(1).total,[0 1]),colorbar,set(gca,'YDir','normal')

     %%%% Figure 8-M: population rate code A1 
    imagesc(Classifier(1).RatePop(1).total,[0 .5]),colorbar,set(gca,'YDir','normal')

    
 L=length(Classifier(2).RastersFm2);
Classifier(2).RatePop(1).total=0;
for s=1:50
    Classifier(2).RatePop(1).total=Classifier(2).RatePop(1).total+Classifier(2).RatePop(1).Mean;
end
    Classifier(2).RatePop(1).total=Classifier(2).RatePop(1).total/L;
    
         %%%% Figure 8-N: population rate code VAF 
    imagesc(Classifier(2).RatePop(1).total,[0 0.5]),colorbar,set(gca,'YDir','normal')
 
    L=length(Classifier(2).RastersFm2);
Classifier(2).RatePop(1).total2=0;
for s=1:50
    Classifier(2).RatePop(1).total2=Classifier(2).RatePop(1).total2+Classifier(2).RatePop(1).Mean2;
end
    Classifier(2).RatePop(1).total2=Classifier(2).RatePop(1).total2/L;

    imagesc(Classifier(2).RatePop(1).total2,[0 0.5]),colorbar,set(gca,'YDir','normal')

       
 L=50;
Classifier(3).RatePop(1).total=0;
for s=1:L
    Classifier(3).RatePop(1).total=Classifier(3).RatePop(1).total+Classifier(3).RatePop(s).Mean;
end
    Classifier(3).RatePop(1).total=Classifier(3).RatePop(1).total/L;
           %%%% Figure 8-O: population rate code SRAF 
imagesc(Classifier(3).RatePop(1).total,[0 .5]),colorbar,set(gca,'YDir','normal')
    
    
  L=150;
Classifier(1).TotalPopulationRate(1).total=0;
for s=1:L
    Classifier(1).TotalPopulationRate(1).total=Classifier(1).TotalPopulationRate(1).total+Classifier(1).TotalPopulationRate(s).Mean;
end
    Classifier(1).TotalPopulationRate(1).total=Classifier(1).TotalPopulationRate(1).total/L;
    imagesc(Classifier(1).TotalPopulationRate(1).total,[0 .5]),colorbar,set(gca,'YDir','normal')
    
    
Fc=[2,4,6,8,11,16,23,32,45,64];
 for k=1:50
A1diagmeanRFm2(k,:)=diag(Classifier(1).RatePop(k).Mean);
  end
A1diagbootR=bootstrp(1000,'mean',A1diagmeanRFm2);
errorbar(Fc,mean(A1diagmeanRFm2).*100,std(A1diagbootR).*100,'b')
xlabel('Fc (Hz)')
ylabel('% correct')
title('A1 population % correct Rate code')
ylim([0 50])


 for k=1:50
VAFdiagmeanRFm2(k,:)=diag(Classifier(2).RatePop(k).Mean);
  end
VAFdiagbootRFm2=bootstrp(1000,'mean',VAFdiagmeanRFm2);
errorbar(Fc,mean(VAFdiagmeanRFm2).*100,std(VAFdiagbootRFm2).*100,'g')
xlabel('Fc (Hz)')
  ylabel('% correct')
  title('VAF population % correct Rate code')
  ylim([0 30])


 for k=1:50
SRAFdiagmeanRFm2(k,:)=diag(Classifier(3).RatePop(k).Mean);
  end
SRAFdiagbootRFm2=bootstrp(1000,'mean',SRAFdiagmeanRFm2);
errorbar(Fc,mean(SRAFdiagmeanRFm2).*100,std(SRAFdiagbootRFm2).*100,'r')
xlabel('Fc (Hz)')
  ylabel('% correct')
  title('SRAF population % correct')
  ylim([0 30])
  
  %%%% Figure 8-L (part2): population rate code A1,VAF,SRAF diagnoal
  errorbar(Fc,mean(A1diagmeanRFm2).*100,std(A1diagbootR).*100,'b','LineWidth', 2)
  hold on 
  errorbar(Fc,mean(VAFdiagmeanRFm2).*100,std(VAFdiagbootRFm2).*100,'g','LineWidth', 2)
  errorbar(Fc,mean(SRAFdiagmeanRFm2).*100,std(SRAFdiagbootRFm2).*100,'r','LineWidth', 2)
%   xlabel('Sound shape')
  xlabel('Fc (Hz)')
  ylabel('% correct')
  title('A1, VAF, SRAF population % correct Rate Code')
  ylim([0 50])
  legend('A1', 'VAF', 'SRAF')
  set(gca,'xscale', 'log')
  
  
   %%% plotting 50 A1 neurons vs %correct with error shade 
 for k=1:50
A1diagRFm2(k,:)=diag(Classifier(1).RatePop(k).Mean);
A1bootRFm2(k,:)=bootstrp(1000,'mean',A1diagRFm2(k,:));
A1eRFm2(k)=std(A1bootRFm2(k,:));
 end
x=1:50;
x=x';
A1yRFm2=mean(A1diagRFm2,2)*100;
A1eRFm2=A1eRFm2'*100;
shadedErrorBar(x,A1yRFm2,A1eRFm2,'b',.5)

[A1ysortRFm2,IsortRFm2]=sort(A1yRFm2,1);
for h=1:50
A1esortRFm2(h)=A1eRFm2(IsortRFm2(h));
end
shadedErrorBar(x,A1ysortRFm2,A1esortRFm2,'b',.5)

title('A1 Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 100])

%%% plotting 50 VAF neurons vs %correct with error shade 
 for k=1:50
VAFdiagRFm2(k,:)=diag(Classifier(2).RatePop(k).Mean);
VAFbootRFm2(k,:)=bootstrp(1000,'mean',VAFdiagRFm2(k,:));
VAFeRFm2(k)=std(VAFbootRFm2(k,:));
 end
x=1:50;
x=x';
VAFyRFm2=mean(VAFdiagRFm2,2)*100;
VAFeRFm2=VAFeRFm2'*100;
shadedErrorBar(x,VAFyRFm2,VAFeRFm2,'g',.5)

[VAFysortRFm2,IsortRFm2]=sort(VAFyRFm2,1);
for h=1:50
VAFesortRFm2(h)=VAFeRFm2(IsortRFm2(h));
end
shadedErrorBar(x,VAFysortRFm2,VAFesortRFm2,'g',.5)

title('VAF Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 30])


%%% plotting 50 SRAF neurons vs %correct with error shade 
 for k=1:50
SRAFdiagRFm2(k,:)=diag(Classifier(3).RatePop(k).Mean);
SRAFbootRFm2(k,:)=bootstrp(1000,'mean',SRAFdiagRFm2(k,:));
SRAFeRFm2(k)=std(SRAFbootRFm2(k,:));
 end
x=1:50;
x=x';
SRAFyRFm2=mean(SRAFdiagRFm2,2)*100;
SRAFeRFm2=SRAFeRFm2'*100;
shadedErrorBar(x,SRAFyRFm2,SRAFeRFm2,'r',.5)

[SRAFysortRFm2,IsortRFm2]=sort(SRAFyRFm2,1);
for h=1:50
SRAFesortRFm2(h)=SRAFeRFm2(IsortRFm2(h));
end
shadedErrorBar(x,SRAFysortRFm2,SRAFesortRFm2,'r',.5)

title('SRAF Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 30])


shadedErrorBar(x,A1yRFm2,A1eRFm2,'b',.5)
hold on 
shadedErrorBar(x,VAFyRFm2,VAFeRFm2,'g',.5)
shadedErrorBar(x,SRAFyRFm2,SRAFeRFm2,'r',.5)
title('A1,VAF,SRAF Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 30])

%%%% Figure 8-P (part2): population rate # of neurons
shadedErrorBar(x,A1ysortRFm2,A1esortRFm2,'b')
hold on 
shadedErrorBar(x,VAFysortRFm2,VAFesortRFm2,'g')
shadedErrorBar(x,SRAFysortRFm2,SRAFesortRFm2,'r')
title('A1,VAF,SRAF Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 50])
    







 for k=1:150
AlldiagRFm2(k,:)=diag(Classifier(1).TotalPopulationRate(k).Mean);
AllbootRFm2(k,:)=bootstrp(1000,'mean',AlldiagRFm2(k,:));
AlleRFm2(k)=std(AllbootRFm2(k,:));
 end
x=1:150;
x=x';
AllyRFm2=mean(AlldiagRFm2,2)*100;
AlleRFm2=AlleRFm2'*100;
shadedErrorBar(x,AllyRFm2,AlleRFm2,'r',.5)

[AllysortRFm2,IsortRFm2]=sort(AllyRFm2,1);
for h=1:150
AllesortRFm2(h)=AlleRFm2(IsortRFm2(h));
end
shadedErrorBar(x,AllysortRFm2,AllesortRFm2,'r')
title('A1,VAF,SRAF Rate Code')
xlabel('Number of Neurons')
ylabel('Percent Correct')
ylim([0 30])


  Fc=[2,4,6,8,11,16,23,32,45,64];
 for k=1:150
AlldiagmeanRFm2(k,:)=diag(Classifier(1).TotalPopulationRate(k).Mean);
  end
AlldiagbootRFm2=bootstrp(1000,'mean',AlldiagmeanRFm2);
errorbar(Fc,(mean(AlldiagmeanRFm2)).*100,(std(AlldiagbootRFm2)).*100,'b*-','LineWidth', 2)
xlabel('Sound Shape Fc (Hz)')
  ylabel('% correct')
  title('All population % correct Rate')
  ylim([0 100])
    set(gca,'xscale', 'log')
    
  
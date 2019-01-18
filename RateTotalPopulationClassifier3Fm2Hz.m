function [Results] = RateTotalPopulationClassifier3Fm2Hz(Classifier,NB,Fm,units)

% NU=length(Classifier(1).TotalRastersFm2);  %Total number of units for specified region
% N=length(Classifier(1).TotalRastersFm2(1).Units);%Number of sounds
% L=length(Classifier(1).TotalRastersFm2(1).Units(1).Trials);%Number of trials 
% Fs=Classifier(1).TotalRastersFm2(1).Units(1).Trials(1).Fs;   %spike train sampling rate
NU=length(Classifier(1).TotalRastersFm2fifty);  %Total number of units for specified region
N=length(Classifier(1).TotalRastersFm2fifty(1).Units);%Number of sounds
L=length(Classifier(1).TotalRastersFm2fifty(1).Units(1).Trials);%Number of trials 
Fs=Classifier(1).TotalRastersFm2fifty(1).Units(1).Trials(1).Fs;   %spike train sampling rate
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
%                 [RAS]=rasterexpand(Classifier(1).TotalRastersFm2(w).Units(n).Trials(imod(n,:,w)),Fm,4);
                [RAS]=rasterexpand(Classifier(1).TotalRastersFm2fifty(w).Units(n).Trials(imod(n,:,w)),Fm,4);
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
%                     [RASv]=rasterexpand(Classifier(1).TotalRastersFm2(w).Units(m).Trials(ival(q,w)),Fm,4);       %Validatoin response - single trial with multiple cycles
                    [RASv]=rasterexpand(Classifier(1).TotalRastersFm2fifty(w).Units(m).Trials(ival(q,w)),Fm,4);       %Validatoin response - single trial with multiple cycles
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
% Last massive change: 06.04.2017
%%% Main %%%

function [HTotCorr,HNoiseCorr]=FirstCorrection(Fracs,Nbins,Y,Words,Nsmpl)

% Preallocating variables
htottemp    = zeros(Nsmpl,1);
hnoisetemp  = zeros(Nsmpl,1);
htot        = zeros(length(Words),length(Fracs));
hnoise      = zeros(length(Words),length(Fracs));
HTotCorr    = zeros(length(Words),1);
HNoiseCorr  = zeros(length(Words),1);


% Select data fractions and calculate entropy
for i = Words
    fprintf(['\n \t calculations for words of length ' num2str(i) '...']);
    for k = 1:1:length(Fracs)
        fprintf(['\n \t\t calculations for spike train of length L/' num2str(Fracs(k)) '...']);
        MaxRangeToPick = Nbins-Nbins/Fracs(k);
        if MaxRangeToPick == 0
            Dataset                 = Y;
            [HtotalRaw, HtotalRate] = Htot(Dataset,i,Ptable(Dataset,i));
            [HnoiseRaw, HnoiseRate] = Hnoise(i,Dataset);
            for j = 1:1:Nsmpl
                htottemp(j,k)       = HtotalRate;
                hnoisetemp(j,k)     = HnoiseRate;
            end
        else
            for j = 1:1:Nsmpl
                StartPoint              = floor(rand(1)*MaxRangeToPick);
                Dataset                 = Y(:,StartPoint:1:StartPoint+floor(Nbins/Fracs(k)));
                [HtotalRaw, HtotalRate] = Htot(Dataset,i,Ptable(Dataset,i));
                [HnoiseRaw, HnoiseRate] = Hnoise(i,Dataset);
                htottemp(j,k)           = HtotalRate;
                hnoisetemp(j,k)         = HnoiseRate;
            end
        end
    end
    htot(i,:)   = mean(htottemp,1);
    hnoise(i,:) = mean(hnoisetemp,1);

    
    % plot data entropies (first correction)
    figure('Name',['1st correction for words of length ' num2str(i)]);
    h.a = axes;
    hold(h.a,'all');
    h.ht = errorbar(Fracs,htot(i,:),std(htottemp,1)/sqrt(Nsmpl),'-sk');
    h.hn = errorbar(Fracs,hnoise(i,:),std(hnoisetemp,1)/sqrt(Nsmpl),'-sm');
    legend('Htot','Hnoise');
    set(h.ht,'MarkerFaceColor','k');
    set(h.hn,'MarkerFaceColor','m');
    axis([0 max(Fracs)+0.5 0.8*min(hnoise(i,:)) 1.2*max(htot(i,:))]);
    xlabel('Data fraction');
    ylabel('Entropy [bits/sec]');
    xticks(Fracs);
    
    % fit
    coeffNames  = {'a0','a1','a2'};
    myfun       = fittype(...
        'a0+a1/x+a2/x^2',...
        'independent','x',...
        'coefficients',coeffNames);
    options     = fitoptions(...
        'method','NonLinearLeastSquares',...
        'StartPoint',[15 1 0.1],...
        'MaxFunEvals',5000,...
        'TolFun',1e-07,...
        'TolX',1e-07,...
        'Lower',[0 -Inf -Inf],...
        'Upper',[+Inf +Inf +Inf]);
    
    %fit for htot
    [ht.cfun,ht.gof]    = fit(1./Fracs',htot(i,:)',myfun,options);
    h.htf               = plot([0 Fracs]',ht.cfun([Inf 1./Fracs]'),'-','LineWidth',1,'Color','r','Parent',h.a);
    
    %fit for hnoise
    [hn.cfun,hn.gof]    = fit(1./Fracs',hnoise(i,:)',myfun,options);
    h.hnf               = plot([0 Fracs]',hn.cfun([Inf 1./Fracs]'),'-','LineWidth',1,'Color','r','Parent',h.a);
    
    % Calculation of the corrected values for the entropies (values for the
    % entropies at 1/fracs=0)
    HTotCorr(i,1)       = ht.cfun(Inf);
    HNoiseCorr(i,1)     = hn.cfun(Inf);


    
    fprintf('\n');
end

clearvars -except HTotCorr HNoiseCorr

end

%%% Subfunctions %%%

% calculate entropy tables
function frequencyCount = Ptable(Data,WordLength)

[Nrep, Nbins] = size(Data);
for k =1:1:Nrep
    frequencyCount{k} = zeros(2.^WordLength,1);
    for j = 1:Nbins-WordLength+1
        word                        = myfasterbin2dec(Data(k,j:(j+WordLength-1)));
        frequencyCount{k}(word+1)   = frequencyCount{k}(word+1)+1;
    end
    frequencyCount{k} = frequencyCount{k}/sum(frequencyCount{k});
end

end

% calculate total entropy
function [HtotalRaw, HtotalRate]= Htot(Data,WordLength,frequencyCount)

[Nrep, Nbins] = size(Data);
for k = 1:1:Nrep
    prob            = frequencyCount{k};
    prob(prob == 0) = [];
    prob(prob == 1) = [];
    if isempty(prob)
        HRaw(k)   = 0;
    else
        HRaw(k)   = -(prob'*log2(prob));
    end
    HRate(k) = HRaw(k)/(WordLength*3e-03);
end
HtotalRaw	= mean(HRaw);
HtotalRate	= mean(HRate);
end

% calculate the entropy tables for noise entropy
function [HnoiseRaw, HnoiseRate] = Hnoise(WordLength, Data)

[Nrep, Nbins] = size(Data);
for j = 1:Nbins-WordLength+1
    frequencyCount = zeros(2^WordLength,1);
    for k = 1:1:Nrep
        word                        = myfasterbin2dec(Data(k,j:(j+WordLength-1)));
        frequencyCount(word+1)      = frequencyCount(word+1)+1;
    end
    prob = frequencyCount/sum(frequencyCount);
    prob(prob == 0) = [];
    prob(prob == 1) = [];
    if isempty(prob)
        HRaw(j)   = 0;
    else
        HRaw(j)   = -(prob'*log2(prob));
    end
    HRate(j)  = HRaw(j)/(WordLength*3e-03);
end
HnoiseRaw   = mean(HRaw);
HnoiseRate  = mean(HRate);
end

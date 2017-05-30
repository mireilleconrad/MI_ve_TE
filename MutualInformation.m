function [MI,Entropy,NoiseEntropy] = MutualInformation(Sim, RandomData, FirstCorr, Nrep, bins, freq, pfail, pspont)
%% Parameters

% Last massive change: 06.04.2017

% Sim         = 0;               % Determine if the results of the simulation are used (1)
% RandomData  = 1;               % Determine if the random generate data are used (1)
% FirstCorr   = 0;               % Determine if the First Correction is applied (1) or not (0)

%bins    = 3;                % ms
%freq    = 20;               % Hz
fracs   = 1:1:5;            % fractions of the dataset to test in the first correction
Nsmpl   = 30;               % #samples for averaging
Words   = [1:10];           % word lengths to compute
%pfail   = 0.7;              % Probability that a spike fails being transmitted
%pspont  = 4*0.07/(1000*3);  % Probability that a spike is spontaneously generated (From non-spike or spike not transmitted)


%% Calculation with the results of the Simulations
if Sim == 1
    
    %Preallocating variables
    htot = zeros(Nsmpl,1);
    hnoise = zeros(Nsmpl,1);

    % what to process
    WTP         = 'NMDA'; 

    % list data files
    dataFiles	= dir(['output/*' WTP '*']);
    nFiles      = length(dataFiles);

    % load parameters
    parameters

    % load conductances
    load('ampa conductance');                               % nS
    load('nmda conductance');                               % nS
    ampa = ampa(deadTime+1:end);
    nmda = nmda(deadTime+1:end);

    % load files and start the processing
    fls = nFiles;
    
    % load data
    fprintf(['\t now dealing with file ' dataFiles(fls).name '\n']);
    load(['output/' dataFiles(fls).name]);
    time    = timeClipped(3:end,:);
    S       = SClipped(3:end,:);
    clear SClipped timeClipped;
    
    % lookup table
    m       = S(:,1);
    h       = S(:,2);
    n       = S(:,3);
    m_iA    = S(:,4);
    h_iA    = S(:,5);
    m_iT    = S(:,6);
    h_iT    = S(:,7);
    Ca      = S(:,8);
    O       = S(:,9);
    V       = S(:,10);
    
    % read gain factor from the filename
    gain(fls) = str2num(dataFiles(fls).name(14:18));
    if gain(fls) == 1
        gainIndex = fls;
    end
    
    % find APs
    [pks,locs]      = findpeaks(S(:,end),'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',3/dt);
    if length(pks) == 0
        frequency(fls) = 0;
    else
        frequency(fls) = 1e03*length(pks)/DTA;
    end
    
    % calculate conductances
    gAMPA   = gain(fls)*1e-06/CS*ampa;                      % mS/cm2
    gNMDA   = gain(fls)*1e-06/CS*nmda;                      % mS/cm2
    gTot    = gAMPA + gNMDA;                                % mS/cm2
    
    % reconstruct sodium currents
    if WTP == 'NMDA'
        iInj    =...
            7/13*gAMPA.*(V-ENa)+...                         % ?A/cm2
            7/13*gNMDA.*(V-ENa).*9.69./(1+0.1688*exp(-0.0717*V));
    elseif WTP == 'CTRL'
        iInj    = 7/13*gTot.*(V-ENa);                       % ?A/cm2
    else
        error('unknow condition');
    end
    ih      = (EK-Eh)/(EK-ENa)*gmax*O.*(V-ENa);             % ih current
    iNa     = gNa*m.^3.*h.*(V-ENa);                         % iNa current
    ET      = 1e03.*(R.*T./(2.*F)).*log(Ca0./Ca);           % calcium reversal potential
    iT      = gT.*m_iT.^M_iT.*h_iT.^N_iT.*(V-ET);           % iT current
    
    % calculate energy load
    energy.EPSC(fls)= -sum(1e-09*iInj*dt);                  % C/cm2
    energy.iNa(fls) = -sum(1e-09*iNa*dt);                   % C/cm2
    energy.ih(fls)  = -sum(1e-09*ih*dt);                    % C/cm2
    energy.iT(fls)  = -sum(1e-09*iT*dt);                    % C/cm2
    
    % binarize and downsample the spike trains
    binarizedSpikeTrain             = zeros(length(S),1);
    binarizedSpikeTrain(locs)       = 1;
    binarizedSpikeTrain             = reshape(binarizedSpikeTrain,length(S)/5,5);    
    binarizedSpikeTrainDownsampled  = zeros((length(binarizedSpikeTrain)-mod(length(binarizedSpikeTrain),3/dt))/(3/dt),5);
    for j = 1:1:5
        for i = 1:1:length(binarizedSpikeTrainDownsampled)
            rng                               = (i-1)*(3/dt)+(1:(3/dt));
            binarizedSpikeTrainDownsampled(i,j)= sum(binarizedSpikeTrain(rng,j));
        end
        binarizedSpikeTrainDownsampled = binarizedSpikeTrainDownsampled > 0;
    end
    
    Y = binarizedSpikeTrainDownsampled';
    [Nrep, Nbins] = size(Y);
end

%% Calculation with random generated Data
if RandomData == 1
    
    [Y,Nbins] = GenerateRandomTransmittedSpikeTrains(bins,freq,pfail,pspont,Nrep);

end

%% Corrections for the Entropy

if FirstCorr == 1
    %First Correction
    fprintf(['\n calculations of the 1st correction ...']);
    [HTotCorr,HNoiseCorr]=FirstCorrection(fracs,Nbins,Y,Words,Nsmpl);
else
    %Entropy Calculation without the first correction
    frequencyCount          = Ptable(Y,Words);
    [HtotalRaw, HtotalRate] = Htot(Y,Words,frequencyCount,bins);
    [HnoiseRaw, HnoiseRate] = Hnoise(Words, Y,bins);
    HTotCorr(:,1)                = HtotalRate;
    HNoiseCorr(:,1)              = HnoiseRate;
end

%Second Correction
% figure('Name','2nd correction');
% fig.a = axes; 
% hold(fig.a,'all');
% fig.htc = plot(1./Words,HTotCorr,'-','LineWidth',2,'Color','r','Parent',fig.a);
% fig.hnc = plot(1./Words,HNoiseCorr,'-','LineWidth',2,'Color','b','Parent',fig.a);
% legend('Htot','Hnoise');
%set(fig.htc,'MarkerFaceColor','r');
%set(fig.hnc,'MarkerFaceColor','b');
% axis([0 max(1./Words)+1 0 1.2*max(HTotCorr)]);
% xticks([0 sort(1./Words)]);
WordsStr = strsplit(rats(sort(1./Words))); WordsStr(length(WordsStr)) = []; WordsStr(1) = [];
% xticklabels([0 WordsStr]);
% xlabel('Inverse of Word lengths');
% ylabel('Entropy [bits/sec]');

% Fit line through total entropy
coeffNames      = {'a','b'};
myfun   = fittype(...
   'a*x+b',...
   'independent','x',...
   'coefficients',coeffNames);
options = fitoptions(...
   'method','NonLinearLeastSquares',...
   'StartPoint',[0.1 22],...
   'MaxFunEvals',5000,...
   'TolFun',1e-07,...
   'TolX',1e-07,...
   'Lower',[-Inf -Inf],...
   'Upper',[+Inf +Inf]);
[ht.fun,gof] = fit(1./Words',HTotCorr,myfun,options);
% fig.htp = plot([0 1./Words],ht.fun([0 1./Words]),'-k','LineWidth',1);

% Fit line through noise entropy
options = fitoptions(...
    'method','NonLinearLeastSquares',...
    'StartPoint',[0.1 1],...
    'MaxFunEvals',5000,...
    'TolFun',1e-07,...
    'TolX',1e-07,...
    'Lower',[-Inf -Inf],...
    'Upper',[+Inf +Inf]);
[hn.fun,gof] = fit(1./Words',HNoiseCorr,myfun,options);
% fig.hnp = plot([0 1./Words],hn.fun([0 1./Words]),'-k','LineWidth',1);

fprintf('\n');

% Mutual Information Calculation
Entropy         = ht.fun(0);
NoiseEntropy    = hn.fun(0);
MI              = Entropy - NoiseEntropy;

end

%% Subfonctions
% calculate entropy tables
function frequencyCount = Ptable(Data,Words)
%
[Nrep, Nbins] = size(Data);
for i = Words
    for k =1:1:Nrep
        frequencyCount{i,k}      = zeros(2^i,1);
        for j = 1:Nbins-i+1
            word                         = myfasterbin2dec(Data(k,j:(j+i-1)));
            frequencyCount{i,k}(word+1)  = frequencyCount{i,k}(word+1)+1;
        end
        frequencyCount{i,k}      = frequencyCount{i,k}/sum(frequencyCount{i,k});
    end
end
end

% calculate total entropy
function [HtotalRaw, HtotalRate]= Htot(Data,Words,frequencyCount,bins)
%
[Nrep, Nbins] = size(Data);
fprintf(['\n calculations for total entropy ...']);
for i = Words
    for k = 1:1:Nrep
        prob            = frequencyCount{i,k};
        prob(prob == 0) = [];
        prob(prob == 1) = [];
        if isempty(prob)
            HRaw(k)   = 0;
            HRate(k)  = 0;
        else
            HRaw(k)   = -(prob'*log2(prob));
            HRate(k)  = HRaw(k)/(i*bins/1000);
        end
    end
HtotalRaw(i)	= mean(HRaw);
HtotalRate(i)	= mean(HRate);
end
end

% calculate the entropy tables for noise entropy
function [HnoiseRaw, HnoiseRate] = Hnoise(Words, Data, bins)
[Nrep, Nbins] = size(Data);
 fprintf(['\n calculations for noise entropy ...']);
for i = Words
    fprintf(['\n \t calculations for words of length ' num2str(i) '...']);
    for j = 1:Nbins-i+1
        frequencyCount = zeros(2^i,1);
        for k = 1:1:Nrep
                word                        = myfasterbin2dec(Data(k,j:(j+i-1)));
                frequencyCount(word+1)      = frequencyCount(word+1)+1;
        end
        prob = frequencyCount/sum(frequencyCount);            
        prob(prob == 0) = [];
        prob(prob == 1) = [];
        if isempty(prob)
            HRaw(j)   = 0;
            HRate(j)  = 0;
        else
            HRaw(j)   = -(prob'*log2(prob));
            HRate(j)  = HRaw(j)/(i*bins/1000);
        end
    end
    HnoiseRaw(i)   = mean(HRaw);
    HnoiseRate(i)  = mean(HRate);                
end
end

% displayed at start
fprintf('now calculating the efficiency curve...\n\n');

% what to process
WTP         = 'NMDA'; 

% list data files
dataFiles	= dir(['output/*' WTP '*']);
nFiles      = length(dataFiles);

% extra parameters
Words       = [1:14];
MaxWordFit  = 13;

% load parameters
parameters

% load conductances
load('ampa conductance');                               % nS
load('nmda conductance');                               % nS
ampa = ampa(deadTime+1:end);
nmda = nmda(deadTime+1:end);

% load files and start the processing
for fls = 1:1:nFiles
    
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
            7/13*gAMPA.*(V-ENa)+...                         % µA/cm2
            7/13*gNMDA.*(V-ENa).*9.69./(1+0.1688*exp(-0.0717*V));
    elseif WTP == 'CTRL'
        iInj    = 7/13*gTot.*(V-ENa);                       % µA/cm2
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
        
    % calculate the entropy tables for total entropy
    fprintf('\n\t\t preprocessing done, calculating probability distributions...\n');
    for i = Words
        fprintf(['\t\t\t calculating for words of length ' num2str(i,'%02.0f') '...']);
        for k = 1:1:5
            frequencyCount{i,k} = zeros(2^i,1);
            for j = 1:length(binarizedSpikeTrainDownsampled)-i+1
                word                        = myfasterbin2dec(binarizedSpikeTrainDownsampled(j:(j+i-1),k));
                frequencyCount{i,k}(word+1) = frequencyCount{i,k}(word+1)+1;
            end
            frequencyCount{i,k} = frequencyCount{i,k}/sum(frequencyCount{i,k});
        end
        fprintf(' done\n');
    end
    
    % calculate total entropy
    for i = Words       
        for k = 1:1:5            
            prob = frequencyCount{i,k};
            prob(prob == 0) = [];
            prob(prob == 1) = [];
            if isempty(prob)                
                HRaw(k)   = 0;
                HRate(k)  = 0;              
            else
                HRaw(k)   = -(prob'*log2(prob));
                HRate(k)  = HRaw(k)/(i*3e-03);
            end
        end        
        Htotal.raw(i)   = mean(HRaw);
        Htotal.rate(i)  = mean(HRate);
        clear HRaw HRate;
    end
    
    % calculate the entropy tables for noise entropy
    fprintf('\n\t\t preprocessing done, calculating probability distributions (noise entropy)...\n');
    for i = Words
        fprintf(['\t\t\t calculating for words of length ' num2str(i,'%02.0f') '...']);
        for j = 1:length(binarizedSpikeTrainDownsampled)-i+1
            frequencyCount = zeros(2^i,1);
            for k = 1:1:5
                    word                        = myfasterbin2dec(binarizedSpikeTrainDownsampled(j:(j+i-1),k));
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
                HRate(j)  = HRaw(j)/(i*3e-03);
            end
        end
        Hnoise.raw(i)   = mean(HRaw);
        Hnoise.rate(i)  = mean(HRate);                
        fprintf(' done\n');
        clear HRaw HRate;
    end
       
    % plot total entropy
    handle.f(fls) = figure('Name',['gain = ' dataFiles(fls).name(14:18)]);
    scl = 1./(Words*3e-03);
    handle.sp(fls,1)    = subplot(211); 
    handle.p(fls,1)     = plot(scl,Htotal.raw,'-or'); hold on;
    handle.p(fls,2)     = plot(scl,Hnoise.raw,'-og');
    ylabel('total entropy (bits)');
    vertscl = max(1.2*max(Htotal.raw),1);
    axis([0 max(scl) 0 vertscl]);
    handle.sp(fls,2)    = subplot(212);
    handle.p(fls,3)     = plot(scl,Htotal.rate,'-or'); hold on;
    handle.p(fls,4)     = plot(scl,Hnoise.rate,'-og'); 
    ylabel('total entropy rate (bits/sec)');
    xlabel('1/word length (sec^{-1})');
    vertscl = max(1.2*max(Htotal.rate),1);
    axis([0 max(scl) 0 vertscl]);
    set([handle.p(fls,:)],'LineWidth',1);
    
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
    [cfun,gof] = fit(scl(1:MaxWordFit)',Htotal.rate(1:MaxWordFit)',myfun,options);
    rng = [0;max(scl)];
    plot(rng,cfun(rng),'-b','LineWidth',2);
%     dummy(fls)      = cfun.b;
    
    % Fit line through noise entropy
    options = fitoptions(...
        'method','NonLinearLeastSquares',...
        'StartPoint',[0.1 1],...
        'MaxFunEvals',5000,...
        'TolFun',1e-07,...
        'TolX',1e-07,...
        'Lower',[-Inf -Inf],...
        'Upper',[+Inf +Inf]);
    [cfun,gof] = fit(scl(1:MaxWordFit)',Hnoise.rate(1:MaxWordFit)',myfun,options);
    plot(rng,cfun(rng),'-b','LineWidth',2); refresh;
%     Hfinal(fls)     = dummy(fls)-cfun.b;
    Hfinal(fls)     = Htotal.rate(10)-Hnoise.rate(10);
    
    % clear variables at the end of the loop
    clear time S gAMPA gNMDA gTot pks locs m h n m_iA h_iA m_iT h_iT Ca O...
        P1 OL V iInj iNa ih binarizedSpikeTrain...
        binarizedSpikeTrainDownsampled frequencyCount prob Htotal Hnoise scl;
    fprintf('\n');
    
end

% convert to appropriate units
elementaryCharge    = 1.602176565e-19;                                          % C
energy.EPSC         = energy.EPSC*CS/(3*elementaryCharge*DTA*1e-03);            % ATP/sec
energy.iNa          = energy.iNa*CS/(3*elementaryCharge*DTA*1e-03);             % ATP/sec
energy.ih           = energy.ih*CS/(3*elementaryCharge*DTA*1e-03);              % ATP/sec
energy.iT           = energy.iT*CS/(6*elementaryCharge*DTA*1e-03);              % ATP/sec

% plot efficiency curve
efficacy        = Hfinal./energy.EPSC;
efficacyNorm    = 100*efficacy/efficacy(gainIndex);
efficacyTotal   = Hfinal./(energy.EPSC+energy.ih+energy.iNa+energy.iT);
efficacyTotalNorm= 100*efficacyTotal/efficacyTotal(gainIndex);
figure; 
plot(gain, Hfinal); 
title('total entropy rate versus gain');
figure;
plot(gain, efficacyNorm);
title('metabolic efficiency versus gain (EPSCs)');
figure;
plot(gain, efficacyTotalNorm);
title('metabolic efficiency versus gain (total)');
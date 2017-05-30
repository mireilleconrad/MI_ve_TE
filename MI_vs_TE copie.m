function MI_vs_TE(MIcalc,TEcalc)

% Last Massive change: 06.04.2017

% Begin global timer and timer for MI calculation
now1 = tic();

% Parameters
rep         = [2 5 10 15 20 30 40 50 75 100 150 200 300 400 500]; % Vector containing the number of repetitions
bins        = 3;                    % ms
freq        = 4;                   % Hz
pin         = freq/1000*bins;       % probability per bin
dur         = 128000;               % ms
Nbins       = floor(dur/bins);      % #bins
pfail       = 0.9;                  % Probability that a spike fails being transmitted
pspont      = 0.8*3/1000;      % Probability that a spike is spontaneously generated (From non-spike or spike not transmitted)
Nmean       = 10;                   % Number of calculation for MI with the Strong method (MI is then the mean over all the calculations)
words       = 7;               % word lengths to compute

saveName = [datestr(datetime('now'),'dd-mm-yyyy-HH:MM') '_freq(' num2str(freq) ')_pfail(' num2str(pfail)... 
    ')_pspont(' num2str(pspont) ')_rep(' num2str(min(rep)) '-' num2str(max(rep)) ')_dur(' num2str(dur) ')'];

% Preallocating variables
MItot           = zeros(Nmean,length(rep));
MItheo1tot      = zeros(Nmean,length(rep));
%MItheo2tot      = zeros(Nmean,length(rep));
TEin_out_mean   = zeros(length(words),length(rep));
TEout_in_mean   = zeros(length(words),length(rep));
STDin_out       = zeros(length(words),length(rep));
STDout_in       = zeros(length(words),length(rep));

%% Mutual Information Calculation

if MIcalc == 1

    % Calculation of Theoretical value of MI
    [~, ~, MItheo1, MItheo2] = MITheory(bins,freq,pfail,pspont);
    MItheo1tot(:,:)          = MItheo1;
    %MItheo2tot(:,:)          = MItheo2;

    % Calculation of MI with the Strong method
    for j = [1:Nmean]
    
        fprintf(['\n calculations for the ' num2str(j) ' time ...']); 
        MI      = zeros(1,length(rep));
        x       = 0;

        for i = rep
  
            fprintf(['\n calculations for Nrep = ' num2str(i) ' ...']);
            x = x+1;
            [MI(x),~,~] = MutualInformation(0, 1, 0, i, bins, freq, pfail, pspont);

        end

    MItot(j,:)        = MI;

    end

    MImean      = mean(MItot);
    MItheo1mean = mean(MItheo1tot);
    %MItheo2mean = mean(MItheo2tot);

    % End timer for MI calculation
    MITime = toc(now1);

end

%% TE calculation
if TEcalc == 1

    % Begin timer for TE calculation
    now2 = tic();

    % Calculation of Tranfer Entropy
    fprintf(['\n calculations of Transfer Entropy ... \n']);
    y = 0;
    
    for k = rep
    fprintf(['\n calculations for Nrep = ' num2str(k) ' ...']);
    y = y+1;
        for l = words
        fprintf(['\n \t calculations for words of length ' num2str(l) '...']);
            [TEin_out_mean(l,y), TEout_in_mean(l,y), STDin_out(l,y), STDout_in(l,y)] = TransferEntropy(k,Nbins,bins,freq,pfail,pspont,l);
        end
    end

    % End timer for TE calculation
    TETime = toc(now2);
    
end

%% Code end

%end global timer
TotalTime = toc(now1);

% plot
h.a = axes;
hold(h.a,'all');
if MIcalc == 1
    if TEcalc == 1
        h.mi        = errorbar(rep, MImean, std(MItot));
        h.mitheo1   = plot(rep, MItheo1mean);
        %h.mitheo2   = plot(rep, MItheo2mean);
        for l = words
            h.teio      = errorbar(rep, TEin_out_mean(l,:), STDin_out(l,:));
            h.teoi      = errorbar(rep, TEout_in_mean(l,:), STDout_in(l,:));
        end
        legend('MI', 'MItheo1','TEin->out','TEout->in');
    else
        h.mi        = errorbar(rep, MImean, std(MItot));
        h.mitheo1   = plot(rep, MItheo1mean);
        %h.mitheo2   = plot(rep, MItheo2mean);
        legend('MI', 'MItheo1');
    end
else
    if TEcalc == 1
         for l = words
            h.teio      = errorbar(rep, TEin_out_mean(l,:), STDin_out(l,:));
            h.teoi      = errorbar(rep, TEout_in_mean(l,:), STDout_in(l,:));
        end
        legend('TEin->out','TEout->in');
    end
end
xlabel('Number of repetitions');
ylabel('[bit/sec]');

%Save figure
saveas(gcf, ['MI_vs_TE/' saveName '.fig']);

%Save values
if MIcalc == 1
    if TEcalc ==1
        save(['MI_vs_TE/' saveName '.mat'],'rep','MImean','MItheo1mean','TEin_out_mean','TEout_in_mean',...
        'MITime','TETime','TotalTime')
    else
        save(['MI_vs_TE/' saveName '.mat'],'rep','MImean','MItheo1mean','MITime','TotalTime')
    end
else
    if TEcalc ==1
        save(['MI_vs_TE/' saveName '.mat'],'rep','TEin_out_mean','TEout_in_mean','TETime','TotalTime')
    end
end

fprintf(['\n Done! \n \n']);

end
function [Y,Nbins] = GenerateRandomTransmittedSpikeTrains(bins,freq,pfail,pspont,Nrep)

% Parameters
%bins    = 3;                % ms
%freq    = 20;               % Hz
p       = freq/1000*bins;   % probability per bin
dur     = 128000;           % ms
% Nrep    = 10;               % #repetitions
Nbins   = floor(dur/bins);  % #bins

% Preallocating variables
Y       = zeros(Nrep,Nbins) > 1;

% Generate input and output
X       = rand(1,Nbins) < p;
Ydel    = rand(Nrep,Nbins) < (1-pfail);
Yadd    = rand(Nrep,Nbins) < pspont;
for i = 1:1:Nrep
    Y(i,:) = X.*Ydel(i,:)+Yadd(i,:);
end

% % plot trains
% figure('Name',['prob. dropping spikes = ' num2str(pfail)...
%     ' (' num2str(Nrep) ' repetitions; input at ' num2str(freq) ' Hz)']);
% axes;
% hold(gca,'all');
% h.input = PlotPSTHStyle(X,Nrep+1,[1 0 0],gca);
% for i = 1:1:Nrep
%     h.output(i) = PlotPSTHStyle(Y(i,:),Nrep+1-i,[0 0 1],gca);
% end
% axis([0 length(X) 0 Nrep+1.5]);
% xlabel('Time [bins]');
% ylabel('Input and repetitions');
% yticks(1:1:Nrep+1);
% for i = 1:1:Nrep
%     str{i} = num2str(Nrep+1-i);
% end
% str{Nrep+1} = 'input';
% yticklabels(str);
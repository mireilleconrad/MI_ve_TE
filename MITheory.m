function [Entropytheo, NoiseEntropytheo, MItheo1, MItheo2] = MITheory(bins,freq,pfail,pspont)

% Last massive change: 06.04.2017

pin = freq/1000*bins;            % probability per bin

p1  = pin*(1-pfail);             % Probability that a spike is transmitted 
p2  = pin*pfail*pspont;          % Probability that a spike is not transmitted and a spike is spontaneously generated
p3  = pin*pfail*(1-pspont);      % Probability that a spike is not transmitted and no spike is spontaneously generated
p4  = (1-pin)*pspont;            % Probability that a "non-spike" spontaneously generates a spike
p5  = (1-pin)*(1-pspont);        % Probability that a "non-spike" does'nt spontaneously generates a spike

% r: Input
pr1 = pin;                       % p(r=1)
pr0 = 1-pin;                     % p(r=0)

% s: output
ps1 = p1+p2+p4;                  % p(s=1)
ps0 = p3+p5;                     % p(s=0)

pr1s1 = p1+p2;                   % p(r=1,s=1)
pr1s0 = p3;                      % p(r=1,s=0)
pr0s1 = p4;                      % p(r=0,s=1)
pr0s0 = p5;                      % p(r=0,s=0)

pr1ks1 = pr1s1/ps1;              % p(r=1|s=1)  k stands for knowing
pr1ks0 = pr1s0/ps0;              % p(r=1|s=0)
pr0ks1 = pr0s1/ps1;              % p(r=0|s=1)
pr0ks0 = pr0s0/ps0;              % p(r=0|s=0)


% MI = H+Hnoise = Sum_r p(r)log(p(r)) + Sum_r,s p(s)p(r|s)log(p(r|s))
%               = Sum_r,s p(r,s)log(p(r,s)/(p(r)p(s)))
% from "Theoretical Neuroscience" P. Dayan and L.F. Abbott, Chapter 4

H       = - pr1*log2(pr1) - pr0*log2(pr0);
Hnoise  = ps1*(pr1ks1*log2(pr1ks1)+pr0ks1*log2(pr0ks1)) + ps0*(pr1ks0*log2(pr1ks0)+pr0ks0*log2(pr0ks0));

Entropytheo         = H/(bins/1000);
NoiseEntropytheo    = Hnoise/(bins/1000);

MItheo1             = Entropytheo + NoiseEntropytheo;
MItheo2             = (pr1s1*log2(pr1s1/(pr1*ps1)) + pr1s0*log2(pr1s0/(pr1*ps0)) + pr0s1*log2(pr0s1/(pr0*ps1)) + pr0s0*log2(pr0s0/(pr0*ps0)))/(bins/1000);


end
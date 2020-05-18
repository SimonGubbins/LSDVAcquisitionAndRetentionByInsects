% Matlab script to update estimates for the duration of infection for
% lumpy skin disease virus in subclinical cattle

% Set the number of chains, etc.
nchains=2;
nburnin=5000;
nsamp=20000;
nthin=2;

%==========================================================================
% INFECTIOUS PERIOD
% Blood PCR (min, max)
D1=[23 25;   % Tuppurainen et al. 2005
    23 25;
    10 12;
    17 19];
D2=[1 4];    % Moeller et al. 2019
D3=[15 17];  % Sohier et al. 2019
D4=[1 6;     % present study
    16 Inf;
    6 10;
    10 14];
durBloodPCR=ParEst_DurEIToo([D1; D2; D3; D4],nchains,nsamp,nburnin,nthin);
%==========================================================================

%==========================================================================
% LATENT PERIOD
% Onset of PCR in blood (min, max)
D1=[3 4;    % Tuppurainen et al. 2005
    2 3;
    2 3];
D2=[3 5];   % Moeller et al. 2019
D3=[4 5];   % Sohier et al. 2019
D4=[11 15;  % presentStudy
    3 5;
    7 9;
    5 7];
onsetBloodPCR=ParEst_DurEIToo([D1; D2; D3; D4],nchains,nsamp,nburnin,nthin);
%==========================================================================

% Save the MCMC samples
save('DurEIToo_MCMCSamples','durBloodPCR','onsetBloodPCR')

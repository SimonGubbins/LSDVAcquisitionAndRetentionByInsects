% Matlab script to update estimates for the duration of infection for
% lumpy skin disease virus in clinical cattle

% Set the number of chains, etc.
nchains=2;
nburnin=5000;
nsamp=20000;
nthin=2;

% Load the MCMC samples for the historic data (see Gubbins 2019)
varload=load('Durations_MCMCSamples');

%==========================================================================
% INFECTIOUS PERIOD (DIFFERENT PROXY MEASURES)
% Blood PCR (min, max)
priorPDF={fitdist([varload.durBloodPCR{1}(:,1);
                   varload.durBloodPCR{2}(:,1)],...
                  'kernel','support','positive');
          fitdist([varload.durBloodPCR{1}(:,2);
                   varload.durBloodPCR{2}(:,2)],...
                  'kernel','support','positive')};        
D=[16 Inf;
   16 Inf;
   16 Inf];
durBloodPCR=ParEst_DurEI(D,priorPDF,nchains,nsamp,nburnin,nthin);

% Skin biopsy (lesions)
priorPDF={fitdist([varload.durSkin{1}(:,1);
                   varload.durSkin{2}(:,1)],...
                  'kernel','support','positive');
          fitdist([varload.durSkin{1}(:,2);
                   varload.durSkin{2}(:,2)],...
                  'kernel','support','positive')};        
D=[12 Inf;
   12 Inf;
   12 Inf];
durSkin=ParEst_DurEI(D,priorPDF,nchains,nsamp,nburnin,nthin);
%==========================================================================

%==========================================================================
% LATENT PERIOD (DIFFERENT PROXY MEASURES)
% Onset of PCR in blood (min, max)
priorPDF={fitdist([varload.onsetBloodPCR{1}(:,1);
                   varload.onsetBloodPCR{2}(:,1)],...
                  'kernel','support','positive');
          fitdist([varload.onsetBloodPCR{1}(:,2);
                   varload.onsetBloodPCR{2}(:,2)],...
                  'kernel','support','positive')};        
D=[3 5;
   3 5;
   3 5];
onsetBloodPCR=ParEst_DurEI(D,priorPDF,nchains,nsamp,nburnin,nthin);

% Onset of skin lesions
priorPDF={fitdist([varload.onsetSkinLesions{1}(:,1);
                   varload.onsetSkinLesions{2}(:,1)],...
                  'kernel','support','positive');
          fitdist([varload.onsetSkinLesions{1}(:,2);
                   varload.onsetSkinLesions{2}(:,2)],...
                  'kernel','support','positive')};        
D=[7 9;
   7 9;
   7 9];
onsetSkinLesions=ParEst_DurEI(D,priorPDF,nchains,nsamp,nburnin,nthin);
%==========================================================================

% Save the MCMC samples
save('DurEI_MCMCSamples',...
     'durBloodPCR','durSkin','onsetBloodPCR','onsetSkinLesions')

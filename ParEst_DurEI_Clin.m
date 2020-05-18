function ParSamp=ParEst_DurEI_Clin(D,priorPDF,nchains,nsamp,nburnin,nthin)
%
% ParSamp=ParEst_DurEI_Clin(D,priorPDF,nchains,nsamp,nburnin,nthin)
%
% Matlab function to estimate parameters gamma distributed parameters
% related to LSDV infection cattle using Bayesian MCMC methods.
%
% Inputs:
% D - array of observed lower and upper limits for period for each animal
% priorPDF - PDF for priors
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs:
% ParSamp - cell array of MCMC samples for each chain

%==========================================================================
% SAMPLE PARAMETER SPACE
% Create cell arrays storing the chain
ParSamp=cell(1,nchains);

% For each chain ...
parfor chain=1:nchains

% Set the initial point for sampling
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)
        k0=unifrnd(1,5);
        mu0=unifrnd(5,15);
        par=[k0; mu0];
        prior=log(pdf(priorPDF{1},par(1)))+log(pdf(priorPDF{2},par(2)));
        CurrL=sum(log(gamcdf(D(:,2),par(1),par(2)/par(1))-...
                      gamcdf(D(:,1),par(1),par(2)/par(1))));
    end
    
% Create the array storing the chain
    ParSampC=zeros(nsamp/nthin,size(par,1)+2);
    iter=1;

% Set the standard deviations for the proposal distribution
    sd=[5 5];

% Set the counter for the number of accepted samples
    n_accept=zeros(size(sd));

% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp)])

% Generate the new parameter set
        for j=1:length(par)
            par_new=par;
            par_new(j)=par(j)+normrnd(0,sd(j));

% Compute the prior and log-likelihood
            prior_new=log(pdf(priorPDF{1},par_new(1)))+...
                      log(pdf(priorPDF{2},par_new(2)));
            NewL=sum(log(gamcdf(D(:,2),par_new(1),par_new(2)/par_new(1))-...
                         gamcdf(D(:,1),par_new(1),par_new(2)/par_new(1))));

% Test whether of not to accept the new parameter set
            u=unifrnd(0,1);
            if isfinite(NewL+prior_new) && ...
               u<min(1,exp((NewL+prior_new)-(CurrL+prior)))
                CurrL=NewL;
                prior=prior_new;
                par=par_new;
                n_accept(j)=n_accept(j)+1;
            end
        end

% Save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end

% Every one hundred samples during burn-in, tune the standard deviation for
% the proposal distribution to ensure an acceptance rate of 20-40%
        for j=1:size(sd,1)
            if samp<=nburnin && mod(samp,100)==1 && n_accept(j)/samp<0.2
                sd(j)=sd(j)/2;
            elseif samp<=nburnin && mod(samp,100)==1 && n_accept(j)/samp>0.4
                sd(j)=2*sd(j);
            end
        end
        
    end

% Store the chain output
    ParSamp{chain}=ParSampC;
    
end
%==========================================================================

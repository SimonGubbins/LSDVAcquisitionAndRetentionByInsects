function ParEst_I2B(sFlag,nF,nchains,nsamp,nburnin,nthin)
%
% ParEst_I2B(sFlag,nF,nchains,nsamp,nburnin,nthin)
%
% Matlab function to implement a Bayesian MCMC scheme to estimate the
% probability of transmission from insect to bovine for LSDV in different
% insect species.
%
% This version uses the Chihota data for transmission from insect to
% bovine.
%
% Inputs:
% sFlag - flag indicating species:
%         1 - Aedes aegypti
%         2 - Culex quinquefasciatus
%         3 - Culicoides nubeculosus
%         4 - Stomoxys calcitrans
% nF - maximum number of insects refeeding (only needed for C. nubeculosus
%      or S. calcitrans)
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% PREPARATORY STUFF
% Set the species name
sppName={'AeAeg', 'CxQui', 'CNubec', 'SCalc'};

% Load the Chihota data
D=load(['ChihotaData_' sppName{sFlag} '.txt']);
if sFlag==3 || sFlag==4
    D(:,2)=nF;
end

% Create the informative priors for thevirus inactivation rate and the
% probability of transmission from bovine to insect
% (NOTE: these are generated by running ParEst_B2I)
varload=load('B2IModel6_MCMCSamples');
PS=[varload.ParSamp{1}(:,[sFlag 5]); varload.ParSamp{2}(:,[sFlag 5])];
priorPDF={fitdist(PS(:,1),'kernel','support','positive'); ...
          fitdist(PS(:,2),'kernel','support',[0 1])};

% Set the number of parameters
if sFlag==1 || sFlag==2
    npar=3;
elseif sFlag==3 || sFlag==4
    npar=4+size(D,1);
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);
%==========================================================================

% For each chain ...
parfor chain=1:nchains

%==========================================================================
% INITIALISE THE CHAIN
% Set the initial scaling factor for the proposal distribution
    sf=(2.38.^2)/npar;
    SIG=eye(npar);
    
% Set the counter for the number of accepted samples
    n_accept=0;

% Create the arrays storing the output for the chain
    ParSampC=zeros(nsamp/nthin,npar+2);
    iter=1;

% Generate the initial parameters for the chain, ensuring they generate a
% finite log likelihood and prior
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)

% Sample an initial set of parameters
        g=unifrnd(0,5);
        bHV=unifrnd(0,1);
        bVH=unifrnd(0,1);
        par=[g; bHV; bVH];
        if sFlag==3 || sFlag==4
            a=unifrnd(0,1);
            nRF=binornd(D(:,2),1-exp(-a.*D(:,1)));
            par=[par; a; nRF];
        end

% Compute the log-likelihood
        [CurrL, prior]=Lhood_I2B(par,D,priorPDF,sFlag);

    end
%==========================================================================

%==========================================================================
% UPDATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Accept: ' num2str(100*n_accept/samp,3) '%'])

% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.01*eye(npar);
        else
            SIGp=sf.*(SIG+0.01*eye(npar));
        end

% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';

% Compute the log likelihood and prior for the new parameter set
        [NewL, prior_new]=Lhood_I2B(par_new,D,priorPDF,sFlag);

% Test whether to accept the new parameter set
        u=unifrnd(0,1);
        if isfinite(NewL+prior_new) && ...
           u<min(1,exp((NewL+prior_new)-(CurrL+prior)))

% Update the counter
            n_accept=n_accept+1;

% Update the covariance matrix for the proposal distribution
            if n_accept==1
                pbar=mean([par par_new],2);
                SIG=cov([par'; par_new']);
            elseif samp<=nburnin && n_accept>1
                pbar_prev=pbar;
                pbar=(n_accept./(n_accept+1)).*pbar_prev+...
                     (1./(n_accept+1)).*par_new;
                SIG=((n_accept-1)./n_accept).*SIG+...
                    (1./n_accept).*(n_accept.*(pbar_prev*pbar_prev')-...
                                    (n_accept+1).*(pbar*pbar')+...
                                    (par_new*par_new'));
            end

% Update the chain
            CurrL=NewL;
            prior=prior_new;
            par=par_new;

        end

% Every one hundred samples during burn-in, tune the scaling factor
% for the proposal distribution to ensure an acceptance rate of 20-40%
        if samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp<0.2
            sf=sf/2;
        elseif samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp>0.4
            sf=2*sf;
        end
%==========================================================================

%==========================================================================
% STORE THE OUTPUT
% After burn in, save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end
%==========================================================================

    end

% Store the chain
    ParSamp{chain}=ParSampC;

end

%==========================================================================
% SAVE THE CHAINS
% Set the filename
fName=[sppName{sFlag} '_I2B_MCMCSamples'];
if sFlag==3 || sFlag==4
    fName=[fName '_N=' num2str(nF)];
end

% Save the outputs
save(fName,'ParSamp','nburnin','nsamp')
%==========================================================================

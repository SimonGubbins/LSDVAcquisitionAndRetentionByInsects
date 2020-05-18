function ParEst_B2I(mFlag,seeds,nchains,nsamp,nburnin,nthin)
%
% ParEst_B2I(mFlag,pFlag,seeds,nchains,nsamp,nburnin,nthin)
%
% Matlab function to implement a Bayesian MCMC scheme to estimate
% the virus inactivation rate (g) and the probability of transmission from
% bovine to insect for LSDV in four insect species.
%
% This version:
% - estimates parameters for all four species (Aedes, Culex, Culicoides and
%   Stomoxys) in a single analysis
% - estimates the probability of transmission when insects feed
%   on clinical (b) and subclinical (rr*b) hosts
%
% Inputs:
% mFlag - flag indicating model to fit:
%         1 - g, b common, rr=1
%         2 - g varies, b common, rr=1
%         3 - g common, b varies, rr=1
%         4 - g, b varies, rr=1
%         5 - g, b, rr common
%         6 - g varies, b, rr common
%         7 - g common, b varies, rr common
%         8 - g, b varies, rr common
%         9 - g, b common, rr varies
%         10 - g varies, b common, rr varies
%         11 - g common, b, rr varies
%         12 - g, b, rr varies
% seeds - vector of seeds for the random number generators used for each
%         chain
% nchains - number of chains to run
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% PREPARATORY STUFF
% Load the data-set
D=load('InsectData.txt');

% Extract the data needed for the analysis:
% 1-species (see sppList)
% 2-status of calf on which insect fed (1-clinical, 2-subclinical)
% 3-days post feeding when tested (insect)
% 4-pos/neg
D=D(:,[1 3 5 8]);

% Set the number of parameters
nParG=[1 4+2];
nParB=[1 4+2];
nParR=[0 1 4+2];
npar=nParG(1+mod(mFlag-1,2))+...
     nParB(1+mod(ceil(mFlag/2)-1,2))+...
     nParR(ceil(mFlag/4));

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);
%==========================================================================

% For each chain ...
parfor chain=1:nchains

%==========================================================================
% INITIALISE THE CHAIN
% Initialise the random number generator
    rng(seeds(chain),'twister');

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

% Sample an initial set of inactivation rates
        if 1+mod(mFlag-1,2)==1
            g=unifrnd(0,5);
        elseif 1+mod(mFlag-1,2)==2
            g=unifrnd(0,5,4,1);
            p=gamfit(g);
            p_g=[p(1); p(1).*p(2)];
        end

% Sample an initial set of probabilities of transmission from clinical
% bovine to insect
        if 1+mod(ceil(mFlag/2)-1,2)==1
            bHV=unifrnd(0,1);
        elseif 1+mod(ceil(mFlag/2)-1,2)==2
            bHV=unifrnd(0,1,4,1);
            p_bHV=betafit(bHV)';
        end

% Sample an initial set of relative risks for subclinical donors compared
% to clinical ones
        if ceil(mFlag/4)==1
            rr=[];
        elseif ceil(mFlag/4)==2
            rr=unifrnd(0,1);
        elseif ceil(mFlag/4)==3
            rr=unifrnd(0,1,4,1);
            p=gamfit(rr);
            p_rr=[p(1); p(1).*p(2)];
        end
        
% Merge the parameters
        par=[g; bHV; rr];
        if 1+mod(mFlag-1,2)==2
            par=[par; p_g];
        end
        if 1+mod(ceil(mFlag/2)-1,2)==2
            par=[par; p_bHV];
        end
        if ceil(mFlag/4)==3
            par=[par; p_rr];
        end

% Compute the log-likelihood
        [CurrL, prior]=Lhood_B2I(par,D,mFlag);

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
        [NewL, prior_new]=Lhood_B2I(par_new,D,mFlag);

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
% COMPUTE DIC AND pD
% Compute the deviance for each sample
Dev=[];
PS=[];
for chain=1:nchains
    Dev=[Dev; -2*ParSamp{chain}(:,end)];
    PS=[PS; ParSamp{chain}(:,1:end-2)];
end

% Compute the mean deviance
Dbar=mean(Dev);

% Compute the deviance at the posterior mean for the parameters
Dhat=-2*Lhood_B2I(mean(PS,1)',D,mFlag);

% Compute the DIC
DIC=2*Dbar-Dhat;

% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================

%==========================================================================
% SAVE THE CHAINS
% Set the filename
fName=['B2IModel' num2str(mFlag) '_MCMCSamples'];

% Save the outputs
save(fName,'ParSamp','seeds','nburnin','nsamp','DIC','pD')
%==========================================================================

function ParEstToo_B2I(cmpt,mFlag,nchains,nsamp,nburnin,nthin)
%
% ParEstToo_B2I(cmpt,mFlag,nchains,nsamp,nburnin,nthin)
%
% Matlab function to implement a Bayesian MCMC scheme to estimate
% the probability of transmission from bovine to insect and duration of
% infection for LSDV in different insect species.
%
% This version fits models for all species simultaneously and explores the
% impact of viral titre in the host on the probability of infection.
%
% Inputs:
% cmpt - compartment used as proxy for infectiousness ('Blood' or 'Skin')
% mFlag - flag indicating model to fit for probability of transmission and
%         virus inactivation rate (g):
%         1 - independent of titre, g common
%         2 - intercept and slope common for all species, g common
%         3 - intercept varies, slope common, g common
%         4 - intercept common, slope varies, g common
%         5 - intercept and slope vary, g common
%         6 - independent of titre, g varies
%         7 - intercept and slope common for all species, g varies
%         8 - intercept varies, slope common, g varies
%         9 - intercept common, slope varies, g varies
%         10 - intercept and slope vary, g varies
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
% 2-days post feeding when tested (insect)
% 3-viral titre in blood or skin of calf on which insect fed
% 4-pos/neg
if strcmp(cmpt,'Blood')==1
    D=D(:,[1 5 6 8]);
elseif strcmp(cmpt,'Skin')==1
    D=D(:,[1 5 7 8]);
end

% Set the number of parameters (related to g,a,b)
if mFlag==1
    npar=1+1+0;
elseif mFlag==2
    npar=1+1+1;
elseif mFlag==3 || mFlag==4
    npar=1+6+1;
elseif mFlag==5
    npar=1+6+6;
elseif mFlag==6
    npar=6+1+0;
elseif mFlag==7
    npar=6+1+1;
elseif mFlag==8 || mFlag==9
    npar=6+6+1;
elseif mFlag==10
    npar=6+6+6;
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);
%==========================================================================

% Set the seeds for the random number generator
seeds=1:nchains;

% For each chain ...
parfor chain=1:nchains

% Initialise the random number generator
    rng(seeds(chain),'twister');

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

% Sample the initial virus inactivation rate(s)
        if mFlag<=5
            g=unifrnd(0,5);
        elseif mFlag>=6
            g=unifrnd(0,5,4,1);
            p_g=gamfit(g);
        end

% Sample the initial dose response parameters
        if mFlag==1
            a=unifrnd(-1,1);
            par=[g; a];
        elseif mFlag==2
            a=unifrnd(-1,1);
            b=unifrnd(-1,1);
            par=[g; a; b];
        elseif mFlag==3
            a=unifrnd(-1,1,4,1);
            b=unifrnd(-1,1);
            [ma,va]=normfit(a);
            par=[g; a; b; ma; va];
        elseif mFlag==4
            a=unifrnd(-1,1);
            b=unifrnd(-1,1,4,1);
            [mb,vb]=normfit(b);
            par=[g; a; b; mb; vb];
        elseif mFlag==5
            a=unifrnd(-1,1,4,1);
            b=unifrnd(-1,1,4,1);
            [ma,va]=normfit(a);
            [mb,vb]=normfit(b);
            par=[g; a; b; ma; va; mb; vb];
        elseif mFlag==6
            a=unifrnd(-1,1);
            par=[g; a; p_g'];
        elseif mFlag==7
            a=unifrnd(-1,1);
            b=unifrnd(-1,1);
            par=[g; a; b; p_g'];
        elseif mFlag==8
            a=unifrnd(-1,1,4,1);
            b=unifrnd(-1,1);
            [ma,va]=normfit(a);
            par=[g; a; b; p_g'; ma; va];
        elseif mFlag==9
            a=unifrnd(-1,1);
            b=unifrnd(-1,1,4,1);
            [mb,vb]=normfit(b);
            par=[g; a; b; p_g'; mb; vb];
        elseif mFlag==10
            a=unifrnd(-1,1,4,1);
            b=unifrnd(-1,1,4,1);
            [ma,va]=normfit(a);
            [mb,vb]=normfit(b);
            par=[g; a; b; p_g'; ma; va; mb; vb];
        end

% Compute the log-likelihood
        [CurrL, prior]=LhoodToo_B2I(par,D,mFlag);

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
        [NewL, prior_new]=LhoodToo_B2I(par_new,D,mFlag);

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
Dhat=-2*LhoodToo_B2I(mean(PS,1)',D,mFlag);

% Compute the DIC
DIC=2*Dbar-Dhat;

% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================

%==========================================================================
% SAVE THE CHAINS
% Set the filename
fName=['B2I_' cmpt 'Model' num2str(mFlag) '_MCMCSamples'];

% Save the outputs
save(fName,'ParSamp','nburnin','nsamp','DIC','pD')
%==========================================================================

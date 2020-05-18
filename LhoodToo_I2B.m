function [logL, prior]=LhoodToo_I2B(par,D,chihotaD,priorPDF)
%
% [logL, prior]=LhoodToo_I2B(par,D,chihotaD,priorPDF)
%
% Matlab function to compute the log likelihood and prior for the
% probability of transmission from insect to bovine for LSDV by Stomoxys
% calcitrans using the Chihota and Sciensano data
%
% Inputs:
% par - vector of model parameters
% D - cell array containing the Sciensano data on transmission to cattle:
%     1 - vector of IDs for recipient cattle
%     2 - vector of times of onset of viraemia (days post initial feed)
%     3 - array of challenge data, columns are:
%         1 - recipient ID
%         2 - time batch first fed on recipient
%         3 - time batch last fed on recipient
%         4 - no. flies in batch
%         5 - status of donor (0-subclinical, 1-clinical)
%     4 - time between feeding on donor and feeding on recipient
% chihotaD - array containing the Chihota data on transmission to cattle,
%            columns are:
%            1-time refeeding post initial feed; 2-no. feeding;
%            3-transmission occurred? (0=no, 1=yes)
% priorPDF - cell array of PDFs for the priors
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior for the parameters

%==========================================================================
% PREPATORY STUFF
% Extract the mechanical transmission parameters
g=par(1);
bHV=par(2);
rr=par(3);
bVH=par(4);

% Extract parameters for the time to onset of viraemia
kE=par(5);     % shape parameter
muE=par(6);    % mean

% Extract the probability of refeeding
a=par(7);
nRF1=floor(par(8:size(D{3},1)+7));
nRF2=floor(par(size(D{3},1)+8:end));
%==========================================================================

%==========================================================================
% PRIORS
% Compute the prior for the mechanical transmission parameters
prior=log(pdf(priorPDF{1},g))+...
      log(pdf(priorPDF{2},bHV))+...
      log(pdf(priorPDF{3},rr))+...
      log(unifpdf(bVH,0,1));

% Add the contribution of the viraemia parameters
prior=prior+log(pdf(priorPDF{4},kE))+...
            log(pdf(priorPDF{5},muE));

% Add the contribution of the probability of refeeding
prior=prior+sum(log(binopdf(nRF1,D{3}(:,4),1-exp(-a.*D{4}))))+...
            sum(log(binopdf(nRF2,chihotaD(:,2),1-exp(-a.*chihotaD(:,1)))))+...
            log(gampdf(a,2.2706,2.2578/2.2706));
%==========================================================================

%==========================================================================
% LIKELIHOOD FOR THE SCIENSANO DATA
% Extract the data
recipList=D{1};
tOnset=D{2};
chalD=D{3};
delta=D{4};

% Initialise the log likelihood
logL=0;

% For each recipient ...
for r=1:length(recipList)

% Extract the challenge data for it
    chalDToo=chalD(chalD(:,1)==recipList(r),:);
    nRFToo=nRF1(chalD(:,1)==recipList(r),:);

% Set the challenge times to days first feeding
    t0=min(chalDToo(:,2));
    chalDToo(:,2)=chalDToo(:,2)-t0;
    chalDToo(:,3)=chalDToo(:,3)-t0;

% Create a vector to store the probability of the outcome for each
% challenge
    p=zeros(size(chalDToo,1),1);
    
% For each challenge ...
    for j=1:size(chalDToo,1)
    
% Set the probability of transmission from bovine to insect, depending on
% the status of the donor
        if chalDToo(j,5)==0
            bHVToo=rr.*bHV;
        elseif chalDToo(j,5)==1
            bHVToo=bHV;
        end    

% Compute the probability of transmission at each feeding in the challenge
        t=chalDToo(j,2):chalDToo(j,3);
        q=zeros(length(t),1);
        for k=1:length(t)
            q(k)=1-(1-bVH.*sum(bHVToo.*exp(-g.*(t(1:k)-t(1)+delta)))).^nRFToo(j);
        end

% If the recipient did not become infected ...
        if isnan(tOnset(r))
            p(j)=prod(1-q);

% If the recipient was infected ...
        elseif ~isnan(tOnset(r))
            for k=1:length(q)
                p(j)=p(j)+q(k)*prod(1-q(1:k-1)).*...
                               (gamcdf(tOnset(r)-t(k),kE,muE/kE)-...
                                gamcdf(tOnset(r)-t(k)-1,kE,muE/kE));
            end
        end
    end

% Update the log likelihood
    if isnan(tOnset(r))
        logL=logL+sum(log(p));
    elseif ~isnan(tOnset(r))
        logL=logL+log(1-(1-prod(p)));
    end

end
%==========================================================================

%==========================================================================
% LIKELIHOOD FOR CHIHOTA DATA
% Calculalate the probability of transmission from insect to bovine for
% each challenge
pTrans=1-(1-bVH.*bHV.*exp(-g.*chihotaD(:,1))).^nRF2;

% Calculate the log likelihood
logL=logL+sum(log(binopdf(chihotaD(:,3),1,pTrans)));
%==========================================================================

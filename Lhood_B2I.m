function [logL, prior]=Lhood_B2I(par,D,mFlag)
%
% [logL, prior]=Lhood_B2I(par,D,mFlag)
%
% Matlab function to compute the log likelihood and prior for the
% probability of transmission from bovine to insect and the duration of
% infection for LSDV in insects.
%
% This version estimates the probability of transmission when insects feed
% on clinical and subclinical hosts.
%
% Inputs:
% par - vector of model parameters
% D - data on insect infection and virus survival; columns are:
%     1-species
%     2-status of calf on which insect fed (1-clinical, 2=subclinical)
%     3-days post feeding when tested
%     4-status (0=negative, 1=positive)
% mFlag - flag indicating model being fit
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior for the parameters

%==========================================================================
% EXTRACT THE PARAMETERS
if mFlag==1
    g=par(1)*ones(4,1);
    bHV=par(2)*ones(4,1);
    rr=ones(4,1);
elseif mFlag==2
    g=par(1:4);
    bHV=par(5)*ones(4,1);
    rr=ones(4,1);
    p_g=par(6:7);
elseif mFlag==3
    g=par(1)*ones(4,1);
    bHV=par(2:5);
    rr=ones(4,1);
    p_bHV=par(6:7);
elseif mFlag==4
    g=par(1:4);
    bHV=par(5:8);
    rr=ones(4,1);
    p_g=par(9:10);
    p_bHV=par(11:12);
elseif mFlag==5
    g=par(1)*ones(4,1);
    bHV=par(2)*ones(4,1);
    rr=par(3)*ones(4,1);
elseif mFlag==6
    g=par(1:4);
    bHV=par(5)*ones(4,1);
    rr=par(6)*ones(4,1);
    p_g=par(7:8);
elseif mFlag==7
    g=par(1)*ones(4,1);
    bHV=par(2:5);
    rr=par(6)*ones(4,1);
    p_bHV=par(7:8);
elseif mFlag==8
    g=par(1:4);
    bHV=par(5:8);
    rr=par(9)*ones(4,1);
    p_g=par(10:11);
    p_bHV=par(12:13);
elseif mFlag==9
    g=par(1)*ones(4,1);
    bHV=par(2)*ones(4,1);
    rr=par(3:6);
    p_rr=par(7:8);
elseif mFlag==10
    g=par(1:4);
    bHV=par(5)*ones(4,1);
    rr=par(6:9);
    p_g=par(10:11);
    p_rr=par(12:13);
elseif mFlag==11
    g=par(1)*ones(4,1);
    bHV=par(2:5);
    rr=par(6:9);
    p_bHV=par(10:11);
    p_rr=par(12:13);
elseif mFlag==12
    g=par(1:4);
    bHV=par(5:8);
    rr=par(9:12);
    p_g=par(13:14);
    p_bHV=par(15:16);
    p_rr=par(17:18);
end
%==========================================================================

%==========================================================================
% COMPUTE THE PRIOR
% Initialise the prior
prior=0;

% Add the contribution to the priors for g
if 1+mod(mFlag-1,2)==1
    prior=prior+log(exppdf(g(1),100));
elseif 1+mod(mFlag-1,2)==2
    prior=prior+sum(log(gampdf(g,p_g(1),p_g(2)/p_g(1))))+...
                log(exppdf(p_g(1),1))+...
                log(exppdf(p_g(2),1));
end    

% Add the contribution to the priors for b
if 1+mod(ceil(mFlag/2)-1,2)==1
     prior=prior+log(unifpdf(bHV(1),0,1));
elseif 1+mod(ceil(mFlag/2)-1,2)==2
     prior=prior+sum(log(betapdf(bHV,p_bHV(1),p_bHV(2))))+...
                 log(exppdf(p_bHV(1),1))+...
                 log(exppdf(p_bHV(2),1));
end

% Add the contribution to the priors for rr
if ceil(mFlag/4)==2
    prior=prior+log(exppdf(rr(1),100));
elseif ceil(mFlag/4)==3
    prior=prior+sum(log(gampdf(rr,p_rr(1),p_rr(2))))+...
                log(exppdf(p_rr(1),1))+...
                log(exppdf(p_rr(2),1));
end
%==========================================================================

%==========================================================================
% COMPUTE THE LIKELIHOOD
% Calculate the probability of each insect being viral DNA positive
clin=(D(:,2)==1);
pPos=NaN(size(D,1),1);
pPos(clin)=bHV(D(clin,1)).*exp(-g(D(clin,1)).*D(clin,3));
pPos(~clin)=rr(D(~clin,1)).*bHV(D(~clin,1)).*exp(-g(D(~clin,1)).*D(~clin,3));

% Calculate the log likelihood
logL=sum(log(binopdf(D(:,4),1,pPos)));
%==========================================================================

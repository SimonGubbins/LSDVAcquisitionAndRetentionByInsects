function [logL, prior]=LhoodToo_B2I(par,D,mFlag)
%
% [logL, prior]=LhoodToo_B2I(par,D,mFlag)
%
% Matlab function to compute the log likelihood and prior for the
% probability of transmission from bovine to insect and the duration of
% infection for LSDV in insects
%
% This version fits models for all species simultaneously and explores the
% impact of viral titre in the host on the probability of infection.
%
% Inputs:
% par - vector of model parameters
% D - data on insect infection and virus survival (columns are: 1-species,
%     2-days post feeding when tested; 3-titre in compartment; 4-status
%     (0=negative, 1=positive)
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
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior for the parameters

% Set the number of vector species and donor animals
nV=4;

% Extract the parameters
if mFlag==1
    g=par(1);
    a=par(2);
elseif mFlag==2
    g=par(1);
    a=par(2);
    b=par(3);
elseif mFlag==3
    g=par(1);
    a=par(2:nV+1);
    b=par(nV+2);
    p_a=par(nV+3:nV+4);
elseif mFlag==4
    g=par(1);
    a=par(2);
    b=par(3:nV+2);
    p_b=par(nV+3:nV+4);
elseif mFlag==5
    g=par(1);
    a=par(2:nV+1);
    b=par(nV+2:2*nV+1);
    p_a=par(2*nV+2:2*nV+3);
    p_b=par(2*nV+4:2*nV+5);
elseif mFlag==6
    g=par(1:nV);
    a=par(nV+1);
    p_g=par(nV+2:nV+3);
elseif mFlag==7
    g=par(1:nV);
    a=par(nV+1);
    b=par(nV+2);
    p_g=par(nV+3:nV+4);
elseif mFlag==8
    g=par(1:nV);
    a=par(nV+1:2*nV);
    b=par(2*nV+1);
    p_g=par(2*nV+2:2*nV+3);
    p_a=par(2*nV+4:2*nV+5);
elseif mFlag==9
    g=par(1:nV);
    a=par(nV+1);
    b=par(nV+2:2*nV+1);
    p_g=par(2*nV+2:2*nV+3);
    p_b=par(2*nV+4:2*nV+5);
elseif mFlag==10
    g=par(1:nV);
    a=par(nV+1:2*nV);
    b=par(2*nV+1:3*nV);
    p_g=par(3*nV+1:3*nV+2);
    p_a=par(3*nV+3:3*nV+4);
    p_b=par(3*nV+5:3*nV+6);
end

% Compute the prior
if mFlag==1
    prior=log(exppdf(g,100))+...
          log(normpdf(a,0,10));
elseif mFlag==2
    prior=log(exppdf(g,100))+...
          log(normpdf(a,0,10))+...
          log(normpdf(b,0,10));
elseif mFlag==3
    prior=log(exppdf(g,100))+...
          sum(log(normpdf(a,p_a(1),p_a(2))))+...
          log(normpdf(b,0,10))+...
          log(normpdf(p_a(1),0,10))+...
          log(exppdf(p_a(2),1));
elseif mFlag==4
    prior=log(exppdf(g,100))+...
          log(normpdf(a,0,10))+...
          sum(log(normpdf(b,p_b(1),p_b(2))))+...
          log(normpdf(p_b(1),0,10))+...
          log(exppdf(p_b(2),1));
elseif mFlag==5
    prior=log(exppdf(g,100))+...
          sum(log(normpdf(a,p_a(1),p_a(2))))+...
          sum(log(normpdf(b,p_b(1),p_b(2))))+...
          log(normpdf(p_a(1),0,10))+...
          log(exppdf(p_a(2),1))+...
          log(normpdf(p_b(1),0,10))+...
          log(exppdf(p_b(2),1));
elseif mFlag==6
    prior=sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
          log(normpdf(a,0,10))+...
          log(exppdf(p_g(1),1))+...
          log(exppdf(p_g(2),1));
elseif mFlag==7
    prior=sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
          log(normpdf(a,0,10))+...
          log(normpdf(b,0,10))+...
          log(exppdf(p_g(1),1))+...
          log(exppdf(p_g(2),1));
elseif mFlag==8
    prior=sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
          sum(log(normpdf(a,p_a(1),p_a(2))))+...
          log(normpdf(b,0,10))+...
          log(exppdf(p_g(1),1))+...
          log(exppdf(p_g(2),1))+...
          log(normpdf(p_a(1),0,10))+...
          log(exppdf(p_a(2),1));
elseif mFlag==9
    prior=sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
          log(normpdf(a,0,10))+...
          sum(log(normpdf(b,p_b(1),p_b(2))))+...
          log(exppdf(p_g(1),1))+...
          log(exppdf(p_g(2),1))+...
          log(normpdf(p_b(1),0,10))+...
          log(exppdf(p_b(2),1));
elseif mFlag==10
    prior=sum(log(gampdf(g,p_g(1),p_g(2)./p_g(1))))+...
          sum(log(normpdf(a,p_a(1),p_a(2))))+...
          sum(log(normpdf(b,p_b(1),p_b(2))))+...
          log(exppdf(p_g(1),1))+...
          log(exppdf(p_g(2),1))+...
          log(normpdf(p_a(1),0,10))+...
          log(exppdf(p_a(2),1))+...
          log(normpdf(p_b(1),0,10))+...
          log(exppdf(p_b(2),1));
end

% Calculate the probability of each insect being viral DNA positive
if mFlag==1
    pPos=exp(-g.*D(:,2))./(1+exp(-a));
elseif mFlag==2
    pPos=exp(-g.*D(:,2))./(1+exp(-(a+b.*D(:,3))));
elseif mFlag==3
    pPos=exp(-g.*D(:,2))./(1+exp(-(a(D(:,1))+b.*D(:,3))));
elseif mFlag==4
    pPos=exp(-g.*D(:,2))./(1+exp(-(a+b(D(:,1)).*D(:,3))));
elseif mFlag==5
    pPos=exp(-g.*D(:,2))./(1+exp(-(a(D(:,1))+b(D(:,1)).*D(:,3))));
elseif mFlag==6
    pPos=exp(-g(D(:,1)).*D(:,2))./(1+exp(-a));
elseif mFlag==7
    pPos=exp(-g(D(:,1)).*D(:,2))./(1+exp(-(a+b.*D(:,3))));
elseif mFlag==8
    pPos=exp(-g(D(:,1)).*D(:,2))./(1+exp(-(a(D(:,1))+b.*D(:,3))));
elseif mFlag==9
    pPos=exp(-g(D(:,1)).*D(:,2))./(1+exp(-(a+b(D(:,1)).*D(:,3))));
elseif mFlag==10
    pPos=exp(-g(D(:,1)).*D(:,2))./(1+exp(-(a(D(:,1))+b(D(:,1)).*D(:,3))));
end

% Calculate the log likelihood
logL=sum(log(binopdf(D(~isnan(pPos),4),1,pPos(~isnan(pPos)))));

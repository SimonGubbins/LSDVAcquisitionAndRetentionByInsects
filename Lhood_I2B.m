function [logL, prior]=Lhood_I2B(par,D,priorPDF,sFlag)
%
% [logL, prior]=Lhood_I2B(par,D,priorPDF,sFlag)
%
% Matlab function to compute the log likelihood and prior for the
% probability of transmission from insect to bovine for LSDV in insects
%
% Inputs:
% par - vector of model parameters
% D - data on transmission to cattle, columns are:
%     1-time refeeding post initial feed; 2-no. feeding; 3-transmission
%     occurred? (0=no, 1=yes)
% priorPDF - cell array of PDFs for the priors
% sFlag - flag indicating species (only required for 3 or 4):
%         1 - Aedes aegypti
%         2 - Culex quinquefasciatus
%         3 - Culicoides nubeculosus
%         4 - Stomoxys calcitrans
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior for the parameters

% Extract the mechanical transmission parameters
g=par(1);
bHV=par(2);
bVH=par(3);
if sFlag==3 || sFlag==4
    a=par(4);
    nRF=floor(par(5:end));
end

% Compute the prior for the mechanical transmission parameters
prior=log(pdf(priorPDF{1},g))+...
      log(pdf(priorPDF{2},bHV))+...
      log(unifpdf(bVH,0,1));
if sFlag==3
      prior=prior+sum(log(binopdf(nRF,D(:,2),1-exp(-a.*D(:,1)))))+...
                  log(unifpdf(a,0,0.4));
elseif sFlag==4
      prior=prior+sum(log(binopdf(nRF,D(:,2),1-exp(-a.*D(:,1)))))+...
                  log(gampdf(a,2.2706,2.2578/2.2706));
end

% Calculalate the probability of transmission from insect to bovine for
% each challenge
if sFlag==1 || sFlag==2
    pTrans=1-(1-bVH.*bHV.*exp(-g.*D(:,1))).^D(:,2);
elseif sFlag==3 || sFlag==4
    pTrans=1-(1-bVH.*bHV.*exp(-g.*D(:,1))).^nRF;
end

% Calculate the log likelihood
logL=sum(log(binopdf(D(:,3),1,pTrans)));

function [LogLikeF,Sigma_F,Sigma_phi]=logl_F(Pars, FitPar, Limits, Min_w, w, F_w, sigma_F_hat)
% Calculate the log-likelihood of F, Sigma_F, and Sigma_phi
% Package: +TimeDelay
% Input  : - A vector of parameters to fit:
%            [Tau, Alpha1, Alpha2, gamma]
%          - A vector of flags indicating which parameter to fit, or to use
%            as constant. The vector is for the following parameters:
%            [Tau, Alpha1, Alpha2, gamma]
%            If NaN, then fit the parameter, otherwise, don't fit this
%            parameter and use the value as the parameter value.
%          - A two column matrix of [lower, upper] bounds on each one of
%            the parameters.
%          - Either [Lower] or [Lower Upper] limit on the angular
%            frequencies to use. Default is 0.
%          - Vector of angular frequency (2*pi*f) corresponding to the
%            observations.
%          - F_w: Vector of the total flux as a function of ang. frequency.
%          - sigma_F_hat - error in F_hat.
% Output : - An array of minus Log likelihood of F for the requested
%            Tau,Alpha1,Alpha2,gamma
%          - Sigma_F, column per parameter Tau/Alpha/gamma
%          - Sigma_phi, column per parameter Tau/Alpha/gamma
% Example: 
% [LogLikeF,Sigma_F,Sigma_phi]=TimeDelay.logl_BPL([1 2],[NaN NaN 0 0],[],[],ResLC.w, ResLC.F_w, ResLC.sigma_F_hat)


if isempty(Limits)
    Limits = [1 1000; 1e-5 1000; 0 1000;  1.5 3.5];
end

if isempty(Min_w)
    Min_w = [0 Inf];
end

if numel(Min_w)==1
    Min_w = [Min_w, Inf];
end


% [Tau, Alpha1, Alpha2, gamma]
Ind = 0;
if ~isnan(FitPar(1))
   Tau = FitPar(1);
else
   Ind = Ind + 1;
   Tau = Pars(Ind);
end

if ~isnan(FitPar(2))
   Alpha1 = FitPar(2);
else
   Ind = Ind + 1;
   Alpha1 = Pars(Ind);    
end
if ~isnan(FitPar(3))
   Alpha2 = FitPar(3);
else
   Ind = Ind + 1;
   Alpha2 = Pars(Ind);    
end
if ~isnan(FitPar(4))
   gamma = FitPar(4);
else
   Ind = Ind + 1;
   gamma = Pars(Ind);    
end



if Tau<Limits(1,1) || Tau>Limits(1,2) || ...
       Alpha1<Limits(2,1) || Alpha1>Limits(2,2) || ...
       Alpha2<Limits(3,1) || Alpha2>Limits(3,2) || ...
       gamma<Limits(4,1) || gamma>Limits(4,2)
    LogLikeF  = Inf;
    Sigma_F   = NaN;
    Sigma_phi = NaN;
else


    w_     = w;
    w_(1)  = 1.0;
    % frequency is column vector
    w_     = w_(:);
    % free parameters should be row vector
    Tau    = Tau(:).';
    Alpha1 = Alpha1(:).';
    Alpha2 = Alpha2(:).';
    gamma  = gamma(:).';

    % sigma_phi is a matrix of (frequrncy/columns, parameters/row):
    Sigma_phi = (Alpha1.^2 + Alpha2.^2 + 2.*Alpha1.*Alpha2.*cos(w.*Tau)) .* (abs(w_).^(-gamma));
    %Sigma_phi = Sigma_phi./(numel(w));
    Sigma_F   = Sigma_phi + sigma_F_hat.^2;

    Power_F = (abs(F_w).^2);

    
    Flag = abs(w)>Min_w(1) & abs(w)<Min_w(2);
    LogLikeF = -0.5.*sum(Power_F(Flag) ./ Sigma_F(Flag)) - 0.5.*sum(log(2.*pi.*Sigma_F(Flag)));
    LogLikeF = -LogLikeF;
end

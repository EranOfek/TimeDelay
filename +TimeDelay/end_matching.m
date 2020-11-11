function [F_t,Slope]=end_matching(T,F_t)
% Apply end-matching to a light curve such that first and last point have the same flux.
% Package: +TimeDelay
% Input  : - Vectot of times.
%          - Vector of fluxes.
% Output : - Vector of fluxes after the end-matching operation.
%          - Removed slope.
% Example: [F_t,Slope]=TimeDelay.end_matching(T,F_t)



Slope = (F_t(end) - F_t(1))./(T(end) - T(1));
PolyPar = [Slope 0];
F_t = F_t - polyval(PolyPar,T);
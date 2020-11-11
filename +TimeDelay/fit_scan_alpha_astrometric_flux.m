function [Res]=fit_scan_alpha_astrometric_flux(t,F_t,x_t,sigma_F,sigma_x,varargin)
% Fit the 2 images astrometric-flux time delay model to observations
% Package: +TimeDelay
% Description:
% Input  : - t : vector of times.
%          - F_t: vector of total flux.
%          - x_t: vector of x position.
%          - y_t: vector of y position.
%          - sigma_F: Error in flux.
%          - sigma_x: Error in position.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Solver' - Either @Util.fit.fminsearch_my | @Util.fit.fminunc_my
%                       Default is @Util.fit.fminunc_my
%            'FitPar' - The list of parameters to fit.
%                       [A0, A1, A2, x0, x1, x2, gamma]
%                       If NaN, then will attempt to fit the parameter.
%                       Default is 
%                       [NaN NaN NaN  NaN NaN NaN  NaN NaN NaN  3]
%            'DefPar' - The list of initial guess to use, or the parameter
%                       value if not fitted. Default is
%                       [2 1   0.66  0   1   -1      3]
%            'Limits' - A two column matrix of [lower, upper] bounds on the
%                       parameters. Default is
%                       [0 5;0 2;0 2;  -1 1; -2.1 2.1; -2.1 2.1;  -1 1;    1.5 3.5]
%            'VecInvTau' - A vector of 1/time_delay to attempt fitting.
%                       Default is (1./100:0.5./100:1./10)
%            'Min_w' - minimum w. Default 2.*pi./100.
%            'Verbose' - Default is true.
% Output : - An output structure containing the following fields:
%            .Tau 
%            .LL_H0
%            .LL_H1
%            .BestPar_H0
%            .BestPar_H1
% Example: Res=TimeDelay.fit_scan_alpha_astrometric_flux
%          ResF=TimeDelay.timedelayed_lc;
% Res=TimeDelay.fit_scan_alpha_astrometric_flux(ResF.T, ResF.F_t,ResF.x_t,ResF.y_t,ResF.eps_F_abs,ResF.eps_x_abs);


NPAR2D = 11 -1;
NPAR1D = 8 -1;

InPar = inputParser;
addOptional(InPar,'Solver',@Util.fit.fminunc_my);
addOptional(InPar,'FitPar',[0 1   0.66  0   1   -1      3]);  % [A0, A1, A2, x0, x1, x2, y0, y1, y2, gamma]
addOptional(InPar,'Limits',[0 20;1e-5 10;0 10;  -1 1; -2.1 2.1; -2.1 2.1;   1.5 3.5]); %  without Tau
addOptional(InPar,'VecA1',logspace(log10(0.01),log10(1),40)); %  (0.5:0.05:1.5));
addOptional(InPar,'VecA2dA1',(0.3:0.05:1));

addOptional(InPar,'Tau',14.7); 
addOptional(InPar,'Min_w',2.*pi./100); 
addOptional(InPar,'Verbose',true);
parse(InPar,varargin{:});
InPar = InPar.Results;


% Input arguments: t,F_t,x_t,y_t,sigma_F,sigma_x

N = length(F_t);
t_step = unique(diff(t));
if numel(t_step)>1
    error('Time series must be equally spaced');
end

freqs = TimeDelay.fft_freq(N, t_step);

sigma_F_hat = sigma_F;
sigma_x_hat = sigma_x;
w = 2.*pi*freqs;

F_w = fft(F_t) ./ sqrt(N);
Gx_t = x_t.*F_t;
Gx_w = fft(Gx_t) ./ sqrt(N);
%Gy_t = y_t.*F_t;% In order 
%Gy_w = fft(Gy_t) ./ sqrt(N);

    

%% Main Fitter 

N = length(F_t);

DFT = fft(eye(N), N, 1) ./ sqrt(N);
DFT_dagger = DFT';
LogZ = sum(log(F_t));
Gamma_1_ = ((DFT * diag(F_t.^2)) * DFT_dagger) * sigma_x_hat^2;



Limits = [InPar.Tau, InPar.Tau; InPar.Limits];
FitParsH1 = [InPar.Tau, InPar.FitPar(1:end)];
FitParsH1(3:4) = NaN;

Na1 = numel(InPar.VecA1);
Na2 = numel(InPar.VecA2dA1);
%Res.LL_H0  = nan(Ntau,1);
Res.A1     = InPar.VecA1;
Res.A2dA1  = InPar.VecA2dA1;
Res.LL_xF  = nan(Na1,Na2);
Res.LL_GF  = nan(Na1,Na2);
Res.LL_F   = nan(Na1,Na2);

Na1.*Na2

for Ia1=1:1:Na1
    for Ia2=1:1:Na2
        
        A1 = InPar.VecA1(Ia1);
        A2 = A1.*InPar.VecA2dA1(Ia2);
        
        [LogL_xF,LogLx_GF,LogL_F] = TimeDelay.logl_x_given_F([A1 A2], FitParsH1, Limits, w, Gx_w, t, F_t, F_w,...
                                                     sigma_F_hat, sigma_x_hat, InPar.Min_w, DFT, DFT_dagger, LogZ, Gamma_1_);
                
        Res.LL_xF(Ia1,Ia2) = LogL_xF;
        Res.LL_GF(Ia1,Ia2) = LogLx_GF;
        Res.LL_F(Ia1,Ia2)  = LogL_F;
        
    end
end


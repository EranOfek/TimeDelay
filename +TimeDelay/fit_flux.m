function Res=fit_flux(t,F_t,sigma_F,varargin)
% Fit the 2 images astrometric-flux time delay model to observations
%
% Example: Res=TimeDelay.fit_flux(ResLC.T,ResLC.F_t,ResLC.sigma_F_hat)


InPar = inputParser;
addOptional(InPar,'Solver',@Util.fit.fminunc_my);
addOptional(InPar,'InputFT',false);  % input is t,F_t (false) or w,F_w (true)

addOptional(InPar,'FitPar',[NaN   NaN              3]);  % [A1, A2, gamma]
addOptional(InPar,'DefPar',[1     0.66             3]);  % [A1, A2, gamma]
addOptional(InPar,'Limits',[1e-5 1000;0 1000;       1.5 3.5]); %  without Tau
addOptional(InPar,'VecInvTau',(1./100:0.5./240:1./10)); %1./(20:1:33)); %1./(10.7:1:18.7)); %(2./100:0.5./100:1./10));  % sign has meaning!
addOptional(InPar,'Min_w',2.*pi./100);
addOptional(InPar,'Verbose',false);
parse(InPar,varargin{:});
InPar = InPar.Results;

InPar.FitPar = [InPar.FitPar];
InPar.DefPar = [InPar.DefPar];

options = optimoptions(@fminunc,'Display','off','Diagnostic','off');
AddPars = {options};




VecTau = 1./InPar.VecInvTau; 
Ntau   = numel(VecTau);


if InPar.InputFT
    w   = t(:);
    F_w = F_t;  % assume F_w is already with ortho normalization
    N   = numel(w);
else
    N     = numel(F_t);
    Tstep = unique(diff(t));
    if numel(Tstep)>1 && range(Tstep)>(10000.*eps)
        error('fit_flux works on evenly spaced data - interpolate data');
    end

    freqs = TimeDelay.fft_freq(N, Tstep(1));
    w = 2.*pi*freqs;
    w = w(:);

    F_w = fft(F_t)./sqrt(N);
end



FitParH0  = [0, NaN, 0, InPar.FitPar(end)];
ParH0     = InPar.DefPar(1);
if isnan(InPar.FitPar(end))
    ParH0 = [ParH0, InPar.DefPar(end)];
end
ParH1     = InPar.DefPar(isnan(InPar.FitPar));

Res.Tau    = VecTau(:);
Res.LL_H0  = nan(Ntau,1);
Res.LL_H1  = nan(Ntau,1);
Res.BestPar_H0 = nan(Ntau,numel(ParH0));
Res.BestPar_H1 = nan(Ntau,numel(ParH1));
Res.ExitFlag_H1 = nan(Ntau,1);
Res.ExitFlag_H0 = nan(Ntau,1);


%Res.LL_H0 = TimeDelay.logl_F(ParH0, FitParH0, InPar.Limits, InPar.Min_w, w, F_w, sigma_F);
Limits    = [0 0;InPar.Limits];

% Pars, FitPar, Limits, Min_w, w, F_w, sigma_F_hat
[Res.BestPar_H0,Res.LL_H0,Res.ExitFlag_H0]=InPar.Solver({@TimeDelay.logl_F, ...
                                         FitParH0, Limits, InPar.Min_w, w, F_w, ...
                                         sigma_F},...
                                        ParH0,AddPars{:});
                                    
for Itau=1:1:Ntau
    if InPar.Verbose
        fprintf('Fitting time delay %d of %d   -  Tau=%f\n',Itau,Ntau,VecTau(Itau));
    end
    
    Limits = [VecTau(Itau), VecTau(Itau); InPar.Limits];
    
    FitParH1 = [VecTau(Itau), InPar.FitPar(1:end)];
    
    [Res.BestPar_H1(Itau,:),Res.LL_H1(Itau),Res.ExitFlag_H1(Itau)]=InPar.Solver({@TimeDelay.logl_F, ...
                                         FitParH1, Limits, InPar.Min_w, w, F_w, ...
                                         sigma_F},...
                                        ParH1,AddPars{:});
    
end
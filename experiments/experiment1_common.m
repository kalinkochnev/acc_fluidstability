% ----------------------------
% This function sets up all parameters that are common across experiments.
% ----------------------------
function [fparams] = experiment1_common(fparams)
	% see eq. 19 for values used
	fparams.Ri = 2;
	fparams.Rp = 2; 
	fparams.Pr = 10;
	fparams.omega = 0.5; % in this paper's 2D case f = frequency, not coriolis force
	fparams.f = 0; % I assume the coriolis force is not considered in 2d case?
	fparams.l = 0;
	fparams.tau = 0.01;


	% background shear flow amplitude
	au = sqrt(2*fparams.Pr*(fparams.Rp - 1)/fparams.Ri);
	fparams.Au = @(t) au*cos(fparams.omega *t);
	fparams.Av = @(t) 0;

	% but for the simple case of Au = au cos(wt) eval integral from 0 to inf analytically
	fparams.Bu = @(t) (au/fparams.omega) * sin(fparams.omega * t);
	fparams.Bv = @(t) 0;
end

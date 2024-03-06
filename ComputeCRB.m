function [crb_medina, crb_tau, crb_eta] = ComputeCRB(nCodeSamp, codeSamp, Fs, snrOut, fc, b)



%% Medina et Al. Implementation

V						= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= V == 0;
V						= (((-1).^(abs(V)))./(V.^2))*2;
V(iSel)				= (pi*pi)/3;
W33_					= real(codeSamp'*V*codeSamp);
W33					= W33_*Fs;

w1_					= real(codeSamp'*codeSamp);
w1						= w1_/Fs;
clear iSel V

Nabla					= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= Nabla == 0;
Nabla					= ((-1).^(abs(Nabla)))./Nabla;
Nabla(iSel)			= 0;
w3						= codeSamp'*Nabla*codeSamp;
clear iSel Nabla

phiFunc				= Fs*(W33 - (abs(w3).^2)/w1);
HTau					= 2*phiFunc/real(codeSamp'*codeSamp);

fim					= ((10.^(snrOut(:)/10)))*(2*Fs*Fs)*(W33_/w1_ - real(w3/w1_).^2);
crb_medina.tau			= 1./fim;


D = diag(1:length(codeSamp));

w2 = 1/Fs^2*(codeSamp'*D*codeSamp);
W22 = 1/Fs^3*(codeSamp'*D^2*codeSamp);

fim = ((10.^(snrOut(:)/10)))*Fs*((2*pi*fc)^2)*(W22 - (w2)^2/w1);
crb_medina.b = 1./fim;
%% Implementation assuming Doppler (b) is known and compensated. 

% V parameter
V						= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= V == 0;
V						= (((-1).^(abs(V)))./(V.^2))*2;
V(iSel)				= (pi*pi)/3;

% Nabla parameter 
Nabla					= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= Nabla == 0;
Nabla					= ((-1).^(abs(Nabla)))./Nabla;
Nabla(iSel)			= 0;


% Coefficients considering the narrow-band case
w1 = (1/Fs)*(codeSamp'*codeSamp);
w3 = codeSamp'*Nabla*codeSamp;
W33 = Fs*codeSamp'*V*codeSamp;
varpi = 1/Fs*codeSamp'*D*Nabla*codeSamp;
W22 = 1/Fs^3*(codeSamp'*D^2*codeSamp);
w2 = 1/Fs^2*(codeSamp'*D*codeSamp);


% Compute FIM and CRB
fim	= ((10.^(snrOut(:)/10)))*(2*Fs*Fs)*((2*pi*fc/Fs).^2 + W33_/w1_ - (4*pi*fc/Fs)*imag(w3/w1_) - real(w3/w1_).^2);
crb_tau.tau		= 1./fim;

%% Implementation assuming Doppler and Time-delay are unknown variables.

wc = 2*pi*fc;
% V parameter
V						= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= V == 0;
V						= (((-1).^(abs(V)))./(V.^2))*2;
V(iSel)				= (pi*pi)/3;

% Nabla parameter 
Nabla					= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= Nabla == 0;
Nabla					= ((-1).^(abs(Nabla)))./Nabla;
Nabla(iSel)			= 0;

% Coefficients considering the narrow-band case
w1 = 1/Fs*(codeSamp')*codeSamp;
w2 = 1/Fs^2*(codeSamp'*D*codeSamp);
w3 = codeSamp'*Nabla*codeSamp;
w4 = 1/Fs*(codeSamp'*D*Nabla*codeSamp);
W22 = 1/Fs^3*(codeSamp'*D^2*codeSamp);
W33 = Fs*codeSamp'*V*codeSamp;
% Define Matrices, narrowband impementation
Q = [   [1j*wc*(1-b),       0,      1] ;
        [   0,          1j*wc,      0]];
w = [w1, w2, w3].';

W = [   [w1,    conj(w2),  conj(w3)];
        [w2,    W22,       conj(w4)];
        [w3,    w4,        W33]];

%% Implementation for CAMSAP

% Parameters
SNR = (10.^(snrOut/10));
tmp = (w1*Fs)./(2*SNR); 

% General FIM and CRB
fim11 = real(W33) + w1*(1-b)^2*(2*pi*fc)^2 + imag(w3)*(1-b)*(4*pi*fc) - real(w3)^2/w1;
fim12 = real(w2)*(1-b)*(2*pi*fc)^2 - imag(conj(varpi))*(2*pi*fc) - ((real(w3)*imag(w2))/(w1))*(2*pi*fc);
fim21 = real(conj(w2))*(1-b)*(2*pi*fc)^2 - imag(conj(varpi))*(2*pi*fc) - ((real(w3)*imag(w2))/(w1))*(2*pi*fc);
fim22 = real(W22)*(2*pi*fc)^2 - imag(w2)^2/w1*(2*pi*fc)^2;

% fim   =  Fs.*[fim11 fim12; fim21 fim22];
fim   =  [fim11 fim12; fim21 fim22];

crb_gen = inv(fim);
crb_eta.tau = (1./(2*SNR))*w1*(crb_gen(1,1));
crb_eta.b = (1./(2*SNR))*w1*(crb_gen(2,2));
end
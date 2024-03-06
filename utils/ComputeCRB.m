function [crbb_tau, crbb_dop, crb_tau , crbCode, phiFunc, fim] = ComputeCRB(nCodeSamp, codeSamp, sampFreq, snrOut, fc, b)


wc      = 2*pi*fc;
alphaP  = ((10.^(snrOut(:)/10)));

V						= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= V == 0;
V						= (((-1).^(abs(V)))./(V.^2))*2;
V(iSel)				= (pi*pi)/3;
W33_					= real(codeSamp'*V*codeSamp);
W33					= W33_*sampFreq;

w1_					= real(codeSamp'*codeSamp);
w1				    = w1_/sampFreq;


Nabla					= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= Nabla == 0;
Nabla					= ((-1).^(abs(Nabla)))./Nabla;
Nabla(iSel)			= 0;
w3						= codeSamp'*Nabla*codeSamp;


phiFunc				= sampFreq*(W33 - (abs(w3).^2)/w1);
HTau					= 2*phiFunc/real(codeSamp'*codeSamp);

fim					= ((10.^(snrOut(:)/10)))*(2*sampFreq*sampFreq)*(W33_/w1_ - real(w3/w1_).^2);
crbCode				= 1./fim;

% Doppler 
w2 = 1/sampFreq^2*(codeSamp'*D*codeSamp);
W22 = 1/sampFreq^3*(codeSamp'*D^2*codeSamp);
fim = ((10.^(snrOut(:)/10)))*sampFreq*((2*pi*fc)^2)*(W22 - (w2)^2/w1);
crb_dop				= 1./fim;

fim				   = ((10.^(snrOut(:)/10)))*(2*sampFreq*sampFreq)*( (2*pi*fc/sampFreq).^2 + W33_/w1_ - (4*pi*fc/sampFreq)*imag(w3/w1_) - real(w3/w1_).^2);
crb_tau				= 1./fim;

% Parameters
D = diag([1:length(codeSamp)]);
SNR = (10.^(snrOut(:)/10));
el1 = 1/(2*pi*fc)^2;

% Coefficients
w1 = 1/sampFreq*codeSamp'*codeSamp;
w3 = codeSamp'*Nabla*codeSamp;
W33 = sampFreq*codeSamp'*V*codeSamp;
varpi = 1/sampFreq*codeSamp'*D*Nabla*codeSamp;
W22 = 1/sampFreq^3*(codeSamp'*D^2*codeSamp);
w2 = 1/sampFreq^2*(codeSamp'*D*codeSamp);


% General FIM and CRB
fim11 = real(W33) + w1*(1-b)^2*(2*pi*fc)^2 + imag(w3)*(1-b)*(4*pi*fc) - real(w3)^2/w1;
fim12 = real(w2)*(1-b)*(2*pi*fc)^2 - imag(conj(varpi))*(2*pi*fc) - ((real(w3)*imag(w2))/(w1))*(2*pi*fc);
fim21 = real(conj(w2))*(1-b)*(2*pi*fc)^2 - imag(conj(varpi))*(2*pi*fc) - ((real(w3)*imag(w2))/(w1))*(2*pi*fc);
fim22 = real(W22)*(2*pi*fc)^2 - imag(w2)^2/w1*(2*pi*fc)^2;

% fim   =  sampFreq.*[fim11 fim12; fim21 fim22];
fim   =  [fim11 fim12; fim21 fim22];

crbb_gen = inv(fim);
crbb_tau = (1./(2*SNR))*w1*(crbb_gen(1,1));
crbb_dop = (1./(2*SNR))*w1*(crbb_gen(2,2));

clear iSel Nabla
clear iSel V
end
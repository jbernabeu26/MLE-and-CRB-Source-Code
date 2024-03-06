% Delay estimation: CRBs (for code only and (code + carrier frequency)) and MSE of MLE of delay based on a real positive amplitude model (so-called modified MLE)
% The oversampling of the MLE is performed via sinc (Shannon theorem)
clear all
close all
addpath('./..')
tic

nMC					= 115e3;

if (nMC < 1e3)
   addDisp			= true;	% additional displays (true, false)
   saveData       = false;
else
   addDisp			= false;	% additional displays (true, false)
   saveData       = true;
end

%%% 1a) Signal generation
for choice = [131,132,133,134,135,31,32,33,34,35]

choice

switch choice

%********************************************************************************************************************************************************************
% NonCircular MLE
%********************************************************************************************************************************************************************

case 31
	sampFreq			= 1;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3;            % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)];		  % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape([codeSamp;zeros(nCodeSamp-1,1)]',[],1);   % NonCircular MLE

case 32
	sampFreq			= 2;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3;            % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)];		  % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape([codeSamp;zeros(nCodeSamp-1,1)]',[],1);   % NonCircular MLE

case 33
	sampFreq			= 5;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/2;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 7; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape([codeSamp;zeros(nCodeSamp-1,1)]',[],1);   % NonCircular MLE

case 34
	sampFreq			= 10;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/2;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 11; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape([codeSamp;zeros(nCodeSamp-1,1)]',[],1);   % NonCircular MLE

case 35
	sampFreq			= 20;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/4;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 14; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape([codeSamp;zeros(nCodeSamp-1,1)]',[],1);   % NonCircular MLE


%********************************************************************************************************************************************************************
% Circular MLE
%********************************************************************************************************************************************************************

case 131
	sampFreq			= 1;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3;            % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)];		  % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape(codeSamp,[],1);                           % Circular MLE

case 132
	sampFreq			= 2;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3;            % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)];		  % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape(codeSamp,[],1);                           % Circular MLE

case 133
	sampFreq			= 5;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/2;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 7; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape(codeSamp,[],1);                           % Circular MLE

case 134
	sampFreq			= 10;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/2;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 11; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape(codeSamp,[],1);                           % Circular MLE

case 135
	sampFreq			= 20;
	fc					= 1540;
	nCode				= 1;
   interpDt       = (1/fc)*(10*2 + 1);% the time interval over which the oversampling is performed is [-interpDt/2,interpDt/2]
   nbPtPerTc      = 120e3/4;          % number of parameter value per 1/fc
	selThreshold	= 0.25;				  % (dB20) selection threshold : parameter values for which the (correlation funtion)^2 is >= Max - selThreshold
   snrOut			= [(40:1:56),(56.5:0.5:78)] - 14; % (dB10)
	PRN				= 1;		           % 1 <= PRN <= 32
	%
	code				= reshape(GenerateCACode(PRN),[],1);
	nChips			= numel(code);
	codeSamp			= reshape(ones(sampFreq,1)*code',[],1);
	nCodeSamp		= numel(codeSamp);
	dataSampIn		= reshape(codeSamp,[],1);                           % Circular MLE

end

%***************************************************************************************
% Every time you start MATLAB, the generator resets itself to the same state. 
% Any script or function that calls the random number functions returns the same result whenever you restart.
% Each time you call rng('shuffle'), it reseeds the generator using a different seed based on the current time.
rng('shuffle') 

selThres				= 10^(selThreshold/20);

if (~exist('nCode'))
	nCode				= 1;
end
if (numel(dataSampIn) == numel(codeSamp))
   circularMle    = true;
else
   circularMle    = false;
end
nSimSampIntObs		= numel(dataSampIn);
isRealDataIn		= sum(abs(imag(dataSampIn))) == 0;
nSimSampIntObs		= numel(dataSampIn);
dataSampSignal		= fft(dataSampIn,nSimSampIntObs);
dataSampSignalC   = conj(dataSampSignal);
normCorr				= sqrt(sum(abs(dataSampIn).^2));
if (1)
	codeCorFunc		= fftshift(ifft(dataSampSignal.*dataSampSignalC))/normCorr;	% MLE correlation function (normalized output noise power)
else
	codeCorFunc		= fftshift(ifft(abs(dataSampSignal).^2))/normCorr;				% MLE correlation function (normalized output noise power)
end
[~,iMax]				= max(codeCorFunc);	% <=> max(abs(codeCorFunc).^2)
if (rem(nSimSampIntObs,2) == 0)
   disp(['max( abs( codeCorFunc(2:iMax-1)    - flipud(conj(codeCorFunc(iMax+1:end)))))      = ', num2str(max(abs(codeCorFunc(2:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end))))))])
else
   disp(['max( abs( codeCorFunc(1:iMax-1)    - flipud(conj(codeCorFunc(iMax+1:end)))))      = ', num2str(max(abs(codeCorFunc(1:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end))))))])
end

if (addDisp)
	figure
	if (rem(nSimSampIntObs,2) == 0)
		plot(abs(codeCorFunc(2:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end))))), grid, axis tight
	else
		plot(abs(codeCorFunc(1:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end))))), grid, axis tight
	end
	title('|codeCorFunc(2:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end)))|')
	ylabel('Amplitude (db20)')
	codeCorFuncD	= codeCorFunc/normCorr;			% normalized correlation function
	[~,iMaxD]		= max(codeCorFuncD);		      % <=> max(abs(codeCorFunc).^2)
	xVal				= ((1:nSimSampIntObs) - iMaxD)/sampFreq;
	xVal				= ((1:nSimSampIntObs) - iMaxD);
	figure
	subplot(2,1,1)
	if (isRealDataIn)
		plot(xVal,codeCorFuncD,'-*k','MarkerSize',4), grid, axis tight
	else
		plot(xVal,real(codeCorFuncD),'-*k',xVal,imag(codeCorFuncD),'-or',xVal,abs(codeCorFuncD),'-db','MarkerSize',4), grid, axis tight
		legend('Real Part', 'Imag Part', 'Modulus')
	end
	xlabel('Time (T_s)')
	ylabel('Amplitude (Lin)')
	ylim([-1.1,1.1])
	title('Code Correlation Function')
	subplot(2,1,2)
	iSel				= find(abs(codeCorFuncD) ~= 0);
	iSel				= iSel(1):iSel(end);
	plot(xVal(iSel),20*log10(abs(codeCorFuncD(iSel))),'-*k','MarkerSize',4), grid, axis tight
	xlabel('Time (T_s)')
	ylabel('Amplitude (db20)')
	ylim([-100,5])
	clear iMaxD xVal iSel codeCorFuncD
end

nSampDataIn       = numel(codeCorFunc);
overSampFac			= (fc*nbPtPerTc)/sampFreq;
sampFreqOut			= sampFreq*overSampFac;

if (overSampFac > 1)
   interpDt       = (ceil((interpDt/2)*sampFreqOut)/sampFreqOut);
   tOut           = reshape((-(interpDt*sampFreqOut):1:(interpDt*sampFreqOut))/sampFreqOut,[],1);
	% A the vicinity of the maximum (real-valued positive real-valued) : real(codeCorFunc.*exp((1j*2*pi*fc)*tOut)) ~= Max*cos((2*pi*fc)*tOut))
	iSel				= cos((2*pi*fc)*tOut) >= (1/selThres);
	tOut				= tOut(iSel);
   disp(['numel(tOut) = ',num2str(numel(tOut)),',   round(numel(iSel)/numel(tOut)) = ', num2str(round(numel(iSel)/numel(tOut)))])
   interpDt       = tOut(end) - tOut(1);
   tIn            = ((1:nSampDataIn) - iMax);    % dataIn are sampled at sampFreq
   matShannon     = repmat(tOut*sampFreq,1,nSampDataIn) - repmat(tIn,numel(tOut),1);
   if (circularMle)
      matShannon	= sincp(matShannon,nSampDataIn);
   else
      matShannon	= sinc(matShannon);
   end
   if (rem(nSampDataIn,2) == 0)
      matShannon(:,1)= 0;
   end
   disp(['size(matShannon) = ',num2str(round(numel(matShannon)*8/(1024^3)*100)/100),' GO'])

   codeCorFunc    = matShannon*codeCorFunc;
   [~,iMax]			= max(codeCorFunc);	% <=> max(abs(codeCorFunc).^2)
   it0				= find(tOut == 0);	
	if (iMax == it0)
		disp(['max( abs( codeCorFunc(1:iMax-1)    - flipud(conj(codeCorFunc(iMax+1:end)))))      = ', num2str(max(abs(codeCorFunc(1:iMax-1) - flipud(conj(codeCorFunc(iMax+1:end))))))])
	else
		disp(['iMax ~= it0, [~,iMax] = max(codeCorFunc) and it0 = find(tOut == 0) => iMax - it0 = ',int2str(iMax - it0)])
		disp(['max( abs( codeCorFunc(1:it0-1)    - flipud(conj(codeCorFunc(it0+1:end)))))      = ', num2str(max(abs(codeCorFunc(1:it0-1) - flipud(conj(codeCorFunc(it0+1:end))))))])
		iMax			= it0;
	end

   if (addDisp)
      figure
		subplot(2,1,1)
      plot(tOut,20*log10(max(abs(codeCorFunc)/normCorr,1e-5)),'-*k','MarkerSize',4),grid,axis tight, 
      xlabel('Time (T_s)')
      ylabel('Amplitude (db20)')
      title('Oversampled Code Correlation Function ')
      clear iMaxD iSel dataSampD xVal
   end
else
   overSampFac    = 1;
   sampFreqOut		= sampFreq*overSampFac;
end

tScale				= tOut;
compPhaseDelay		= exp((1j*2*pi*fc)*tScale);
corFunc				= real(codeCorFunc.*compPhaseDelay);
[~,iMax2]			= max(corFunc);				% <=> max((corFunc).^2)
% delayTrue		   = iMax2/sampFreqOut;
delayTrue		   = tOut(iMax2);
if (iMax2 ~= iMax)
	warning(['Max of codeCorFunc and max of corFunc do not coincide : dInd = ',int2str(iMax2-iMax),...
                                                                    ', dt = ',num2str((iMax2-iMax)/sampFreqOut), ' T_s', ...
                                                                    ', dd = ',num2str((iMax2-iMax)/sampFreqOut*300*1000),' mm'])
end

if (addDisp)
	subplot(2,1,2)
	plot(tScale,20*log10(max(abs(codeCorFunc)/normCorr,1e-5)),'-k',tScale,20*log10(max(abs(corFunc)/normCorr,1e-5)),'-*r','MarkerSize',4),grid,axis tight, 
	xlabel('Time (T_s)')
	ylabel('Amplitude (db20)')
	ylim([-9,1])
	title('Oversampled (Code Correlation Function) vs (Correlation Function)')
   clear iMaxD iSel dataSampD xVal
end

drawnow

% keyboard

dataSampSignal    = repmat(dataSampSignal,1,nMC);
dataSampSignalC   = repmat(dataSampSignalC,1,nMC);
compPhaseDelay    = repmat(compPhaseDelay,1,nMC);

amp					= sqrt((10.^(snrOut/10))/real(codeSamp'*codeSamp));
mse					= zeros(numel(amp),1);
delayEst				= zeros(nMC,numel(amp));
maxEst				= zeros(nMC,numel(amp));

%%% Oversampled observation and MLE
for iAmp	= 1:numel(amp)
   disp([int2str(round(iAmp/numel(amp)*100)), ' %'])

   %%% 2a) Noise generation
   dataSampNoise		= fft((randn(nSimSampIntObs,nMC) + 1j* randn(nSimSampIntObs,nMC))/sqrt(2),nSimSampIntObs);      % arr(nSimSampIntObs,nMC)
   sigCorFunc			= fftshift(ifft((dataSampSignal*amp(iAmp) + dataSampNoise).*dataSampSignalC,[],1),1)/normCorr;  % arr(nSimSampIntObs,nMC)
   clear dataSampNoise

   %%% 2b) MLE oversampling
   if (overSampFac > 1)
      sigCorFunc        = matShannon*sigCorFunc;
   end

   sigCorFunc				= real(sigCorFunc.*compPhaseDelay);
   [maxL,iMaxL]			= max(sigCorFunc,[],1);		% <=> max((sigCorFunc).^2) & sigCorFunc > 0
   clear sigCorFunc
   % delayEst(:,iAmp)	= iMaxL(:)/sampFreqOut;
   delayEst(:,iAmp)	   = tScale(iMaxL(:));
   maxEst(:,iAmp)       = maxL(:);
   critSel              = maxEst(:,iAmp) > 0;
	mse(iAmp)				= sum(((delayEst(:,iAmp) - delayTrue).*critSel).^2)/sum(critSel);
end
clear dataSampSignal dataSampSignalC compPhaseDelay matShannon


crb					= mse*NaN;
crbCode				= mse*NaN;
if (saveData)
   if (circularMle)
      fileName    = 'CrbVsMseModMleSincOptCirc_';
   else
      fileName    = 'CrbVsMseModMleSincOptNonCirc_';
   end
   save([fileName,int2str(choice)],'snrOut','fc','nMC','nCode','sampFreq','interpDt','nbPtPerTc','selThreshold','PRN','codeSamp','tScale','codeCorFunc','corFunc','delayTrue','delayEst','maxEst','mse','crb','crbCode')
end

V						= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= V == 0;
V						= (((-1).^(abs(V)))./(V.^2))*2;
V(iSel)				= (pi*pi)/3;
W33_					= real(codeSamp'*V*codeSamp);
W33					= W33_*sampFreq;

w1_					= real(codeSamp'*codeSamp);
w1						= w1_/sampFreq;
clear iSel V

Nabla					= reshape(0:(nCodeSamp-1),[],1)*ones(1,nCodeSamp) - ones(nCodeSamp,1)*(0:(nCodeSamp-1));
iSel					= Nabla == 0;
Nabla					= ((-1).^(abs(Nabla)))./Nabla;
Nabla(iSel)			= 0;
w3						= codeSamp'*Nabla*codeSamp;
clear iSel Nabla

phiFunc				= sampFreq*(W33 - (abs(w3).^2)/w1);
HTau					= 2*phiFunc/real(codeSamp'*codeSamp);

fim					= ((10.^(snrOut(:)/10)))*(2*sampFreq*sampFreq)*(W33_/w1_ - real(w3/w1_).^2);
crbCode				= 1./fim;

fim					= ((10.^(snrOut(:)/10)))*(2*sampFreq*sampFreq)*((2*pi*fc/sampFreq).^2 + W33_/w1_ - (4*pi*fc/sampFreq)*imag(w3/w1_) - real(w3/w1_).^2);
crb					= 1./fim;

if (saveData)
   if (circularMle)
      fileName    = 'CrbVsMseModMleSincOptCirc_';
   else
      fileName    = 'CrbVsMseModMleSincOptNonCirc_';
   end
   save([fileName,int2str(choice)],'snrOut','fc','nMC','nCode','sampFreq','interpDt','nbPtPerTc','selThreshold','PRN','codeSamp','tScale','codeCorFunc','corFunc','delayTrue','delayEst','maxEst','mse','crb','crbCode')
end


figure
plot(snrOut,10*log10(mse),'-*k',snrOut,10*log10(crbCode),'-or',snrOut,10*log10(crb),'-ob'), axis tight , grid
xlabel('SNR OUT (dB10)')
ylabel('MSE (dB10)')
legend(['MLE (F_s = ',int2str(sampFreq),')'],['CRB CODE (F_s = ',int2str(sampFreq),')'],['CRB (F_s = ',int2str(sampFreq),')'])

end % choice

toc



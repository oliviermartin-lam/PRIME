function trs = modelTransferFunctions(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

%% 1\ HO LOOP
%1.1. Define the frequency vectors
trs.holoop.tf.freq = linspace(0,trs.wfs.nExp/2,trs.wfs.nExp/2)*trs.holoop.freq/trs.wfs.nExp;
%1.2 wfs TF
trs.holoop.tf.wfs    = wfsTransferFunction(trs.holoop.tf.freq,trs.holoop.freq);
%1.3. Lag tf
trs.holoop.tf.lag = delayTransferFunction(trs.holoop.tf.freq,trs.holoop.freq,trs.holoop.lat*trs.holoop.freq);
%1.4 servo tf
trs.holoop.tf.servo = servoTransferFunction(trs.holoop.tf.freq,trs.holoop.freq,trs.holoop.tf.num,trs.holoop.tf.den);
%1.5. open loop tf
trs.holoop.tf.ol = trs.holoop.tf.wfs.*trs.holoop.tf.servo.*trs.holoop.tf.lag;
%1.6. closed loop tf
trs.holoop.tf.ctf = trs.holoop.tf.ol./(1+ trs.holoop.tf.ol);
%1.7. rejection tf
trs.holoop.tf.rtf =1 - trs.holoop.tf.ctf;
%1.8. noise transfer function
trs.holoop.tf.ntf = squeeze(trs.holoop.tf.servo./(1+trs.holoop.tf.ol));
trs.holoop.tf.pn  = (trapz(trs.holoop.tf.freq,abs(trs.holoop.tf.ntf).^2)*2/trs.holoop.freq);

%% 2\ TT LOOP
%2.1. Define the frequency vectors
trs.ttloop.tf.freq = linspace(0,trs.tipTilt.nExp/2,trs.tipTilt.nExp/2)*trs.ttloop.freq/trs.tipTilt.nExp;%logspace(-3,log10(0.5*trs.ttloop.freq),trs.tipTilt.nExp/2);
%2.2 wfs TF
trs.ttloop.tf.wfs    = wfsTransferFunction(trs.ttloop.tf.freq,trs.ttloop.freq);
%2.3 TT Lag tf
trs.ttloop.tf.lag = delayTransferFunction(trs.ttloop.tf.freq,trs.ttloop.freq,trs.ttloop.lat*trs.ttloop.freq);
%2.4 TT servo tf
trs.ttloop.tf.servo = servoTransferFunction(trs.ttloop.tf.freq,trs.ttloop.freq,trs.ttloop.tf.num,trs.ttloop.tf.den);
%2.5 open loop tf
trs.ttloop.tf.ol = trs.ttloop.tf.wfs.*trs.ttloop.tf.servo.*trs.ttloop.tf.lag;
%2.6 closed loop tf
trs.ttloop.tf.ctf =trs.ttloop.tf.ol./(1+trs.ttloop.tf.ol);
%2.7 rejection tf
trs.ttloop.tf.rtf =1 - trs.ttloop.tf.ctf;
%2.8 noise transfer function
trs.ttloop.tf.ntf = squeeze(trs.ttloop.tf.servo./(1+trs.ttloop.tf.ol));
trs.ttloop.tf.pn  = (trapz(trs.holoop.tf.freq,abs(trs.ttloop.tf.ntf).^2)*2/trs.ttloop.freq);


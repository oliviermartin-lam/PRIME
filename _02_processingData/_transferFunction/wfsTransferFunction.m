function h = wfsTransferFunction(freq,fSamp)
h    = sinc(freq/fSamp/pi).*exp(-1i*freq/fSamp);
end

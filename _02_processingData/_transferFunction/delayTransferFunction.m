function h = delayTransferFunction(freq,fSamp,delay)
    h    = exp(-2i*freq*delay/fSamp);
end

function h = servoTransferFunction(freq,fSamp,num,den)
z = exp(-2i*pi*freq/fSamp);
ho_num = num(1) + num(2)*z + num(3)*z.^2 + num(4)*z.^3;
ho_den = 1+ den(1)*z + den(2)*z.^2 + den(3)*z.^3;
h = ho_num./ho_den;
end
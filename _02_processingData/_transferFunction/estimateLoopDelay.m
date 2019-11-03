%{
------------HEADER-----------------
Objective          ::  Estimate loops delays
INPUT VARS
wssmprg              :: Header field from the TRS structurue
OUTPUT VARS
r0               :: HO loop delay in seconds
L0               :: TT loop delay in seconds
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [delay,delay_tt] = estimateLoopDelay(wssmprg,aoMode)
inputs = inputParser;
inputs.addRequired('wssmprg',@isnumeric);
inputs.addRequired('aoMode',@ischar);

%Clock-dependent delays from Microgate NGWFC DDR document - verified
%through daytime tests

switch wssmprg
    case 0
        t0 = 8000e-6;
    case 1
        t0 = 2600e-6;
    case 2
        t0 = 850e-6;
    case 3
        t0 = 700e-6;
    case 4
        t0 = 28000e-6;
    case 5
        t0 = 9500e-6;
    case 6
        t0 = 2000e-6;
    case 7
        t0 = 1350e-6;
end

% Total delays
t8 = 40e-6;
delay = t0+t8;

if strcmp(aoMode,'NGS')
    t9 = 1e-6;
    t10 = 2e-6;
    textra_tt = 200e-6;% retrieved from rejection transfer function fit
    delay_tt = t0 + t9 + t10 + textra_tt;
else
    delay_tt = 1e-3;%STRAP is much faster
end

%Vand Dam SPIE 2004: 2.13ms and 1.65ms



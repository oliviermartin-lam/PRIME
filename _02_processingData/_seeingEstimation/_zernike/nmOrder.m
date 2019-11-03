%{
------------HEADER-----------------
Objective          :: Get the radial and azimuthal degree from the j-index.
 
INPUT VARS
jindex              :: Noll Zernike j-index

OUTPUT VARS
n,m            :: Zernike radial and azimuthal degrees
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}


function [n,m] = nmOrder(jindex)
% [n,m] = nmOrder(i,)
%Give the radial and azimutal order of the ith zernike
n = floor((-1.+sqrt(8*(jindex-1)+1))/2);
p = (jindex-(n.*(n+1))/2);
k = mod(n,2);
m = floor((p+k)/2)*2 - k;
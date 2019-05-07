function [ B ] = Bw( P )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
B0=1;
cw=3e-06;
Pref=2800;
B=B0*exp(-cw*(P-Pref));

end


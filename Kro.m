function [ K ] = Kro(Sw)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
Kromax=0.7;
Sor=0.3;
Swr=0.25;
Somin=Sor;
Somax=1-Swr;
So=1-Sw;
if So<=Somin
    K=0;
elseif So>=Somax
    K=Kromax;
else
    K=Kromax*((So-Sor)/(1-Sor-Swr))^3;
end


end


function [ K ] = Krw( Sw )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

Krwmax=0.08;
Sor=0.3;
Swr=0.25;
Swmin=Swr;
Swmax=1-Sor;
if Sw<=Swmin
    K=0;
elseif Sw>=Swmax
    K=Krwmax;
else
    K=Krwmax*((Sw-Swr)/(1-Sor-Swr))^2;
end

end


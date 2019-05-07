%%%%%%%%%%%curve fitting the oil properties%%%%%%%%%%%%%%

P=[400:400:5600];
Bo=[1.012 1.009 1.005 1.001 0.996 0.99 0.988 0.985 0.98 0.975 0.97 0.965 0.96 0.955];
uo=[1.17 1.14 1.11 1.08 1.06 1.03 1 0.98 0.95 0.94 0.92 0.91 0.9 0.89];
B=polyfit(P,Bo,2);
u=polyfit(P,uo,2);
Bo_d=polyval(B,P);
uo_d=polyval(u,P);
plot(P,Bo,'o',P,Bo_d);
figure();
plot(P,uo,'o',P,uo_d);

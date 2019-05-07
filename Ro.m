function [ Res_o ] = Ro( P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,increment_P,increment_Sw,Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


Pbhp=2900;
dx=30;
dy=30;
dz=30;
nx=15;
ny=15;
Qi=300*5.615;
Vp=dx*dy*dz;
uw=1;
e1=10^-6;


if strcmp(increment_P,'+i')==1
    P_v02d(i,j+1)=P_v02d(i,j+1)+e1;
elseif strcmp(increment_P,'-i')==1
    P_v02d(i,j-1)=P_v02d(i,j-1)+e1;
elseif strcmp(increment_P,'+j')==1
    P_v02d(i+1,j)=P_v02d(i+1,j)+e1;
elseif strcmp(increment_P,'-j')==1
    P_v02d(i-1,j)=P_v02d(i-1,j)+e1;
elseif strcmp(increment_P,'ij')==1
    P_v02d(i,j)=P_v02d(i,j)+e1;  
elseif strcmp(increment_P,'none')==1
    P_v02d=P_v02d;
end

if strcmp(increment_Sw,'+i')==1
    Sw_v02d(i,j+1)=Sw_v02d(i,j+1)+e1;
elseif strcmp(increment_Sw,'-i')==1
    Sw_v02d(i,j-1)=Sw_v02d(i,j-1)+e1;
elseif strcmp(increment_Sw,'+j')==1
    Sw_v02d(i+1,j)=Sw_v02d(i+1,j)+e1;
elseif strcmp(increment_Sw,'-j')==1
    Sw_v02d(i-1,j)=Sw_v02d(i-1,j)+e1;
elseif strcmp(increment_Sw,'ij')==1
    Sw_v02d(i,j)=Sw_v02d(i,j)+e1; 
elseif strcmp(increment_Sw,'none')==1
    Sw_v02d=Sw_v02d;
end

%%%%%%%New tranmsissibility calculation%%%
if strcmp(increment_P,'none')~=1 || strcmp(increment_Sw,'none')~=1
[ Tx_o,Ty_o,Tx_w,Ty_w ] = Trans(P_v02d,Sw_v02d,Tx_g,Ty_g);
end
    

%%%%%oil residual
if j==1
    R1_o=Tx_o(i,j)*(0-P_v02d(i,j));
else
R1_o=Tx_o(i,j)*(P_v02d(i,j-1)-P_v02d(i,j));
end

if j==nx
    R2_o=Tx_o(i,j+1)*(0-P_v02d(i,j));
else
R2_o=Tx_o(i,j+1)*(P_v02d(i,j+1)-P_v02d(i,j));
end

if i==1
    R3_o=Ty_o(i,j)*(0-P_v02d(i,j));
else
R3_o=Ty_o(i,j)*(P_v02d(i-1,j)-P_v02d(i,j));
end

if i==ny
    R4_o=Ty_o(i+1,j)*(0-P_v02d(i,j));
else
R4_o=Ty_o(i+1,j)*(P_v02d(i+1,j)-P_v02d(i,j));
end

if i==1&&j==1%% injector well
R5_o=0;
else
    R5_o=0;
end

if i==15&&j==15 %%producer well
   
    WIf_o=Kro(Sw_v02d(i,j))/(Uo(P_v02d(i,j))*Bo(P_v02d(i,j)));
    R6_o=WI_g*WIf_o*(P_v02d(i,j)-Pbhp);
else
    R6_o=0;
end

dphi=(poro(P_v02d(i,j))-poro(P_n2d(i,j)))/(P_v02d(i,j)-P_n2d(i,j));
bo_n1=1/Bo(P_v02d(i,j));
bo_n=1/Bo(P_n2d(i,j));
db_o=(bo_n1-bo_n)/(P_v02d(i,j)-P_n2d(i,j));

R7_o=(Vp/dt)*((-poro(P_v02d(i,j))*bo_n1)*(Sw_v02d(i,j)-Sw_n2d(i,j))+(1-Sw_n2d(i,j))*((dphi*bo_n1+poro(P_n2d(i,j))*db_o))*(P_v02d(i,j)-P_n2d(i,j)));
Res_o=R1_o+R2_o+R3_o+R4_o-R5_o-R6_o-R7_o;




end


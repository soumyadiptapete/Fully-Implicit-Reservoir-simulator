function [ Res_w ] = Rw( P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,increment_P,increment_Sw,Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)
%UNTITLED2 Summary of this function goes here
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
e1=1e-06;

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
    

%%%%%%%%%%%%water residual%%%
if j==1
R1_w=Tx_w(i,j)*(0-P_v02d(i,j));   
else
R1_w=Tx_w(i,j)*(P_v02d(i,j-1)-P_v02d(i,j));
end
if j==nx
    R2_w=Tx_w(i,j+1)*(0-P_v02d(i,j));
else
R2_w=Tx_w(i,j+1)*(P_v02d(i,j+1)-P_v02d(i,j));
end

if i==1
    R3_w=Ty_w(i,j)*(0-P_v02d(i,j));
else
R3_w=Ty_w(i,j)*(P_v02d(i-1,j)-P_v02d(i,j));
end

if i==ny
    R4_w=Ty_w(i+1,j)*(0-P_v02d(i,j));
else
R4_w=Ty_w(i+1,j)*(P_v02d(i+1,j)-P_v02d(i,j));
end

if i==1&&j==1%% injector well
R5_w=-Qi;
else
    R5_w=0;
end

if i==15&&j==15 %%producer well
    
    WIf_w=Krw(Sw_v02d(i,j))/(uw*Bw(P_v02d(i,j)));
    R6_w=WI_g*WIf_w*(P_v02d(i,j)-Pbhp);
else
    R6_w=0;
end


bw_n1=1/Bw(P_v02d(i,j));
bw_n=1/Bw(P_n2d(i,j));
dbw=(bw_n1-bw_n)/(P_v02d(i,j)-P_n2d(i,j));
dphi=(poro(P_v02d(i,j))-poro(P_n2d(i,j)))/(P_v02d(i,j)-P_n2d(i,j));
R7_w=(Vp/dt)*((poro(P_v02d(i,j))*bw_n1)*(Sw_v02d(i,j)-Sw_n2d(i,j))+(Sw_n2d(i,j)*bw_n1*dphi+Sw_n2d(i,j)*poro(P_n2d(i,j))*dbw)*(P_v02d(i,j)-P_n2d(i,j)));


Res_w=R1_w+R2_w+R3_w+R4_w-R5_w-R6_w-R7_w;


end


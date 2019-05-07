function [ Tx_o,Ty_o,Tx_w,Ty_w ] = Trans(P_v02d,Sw_v02d,Tx_g,Ty_g)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ny=15;
nx=15;
%%% X -direction
Tx_o=zeros(ny,nx+1);
Tx_w=zeros(ny,nx+1);
for i=1:ny
    for j=2:nx
        
        %%% potential of flow from i-1 to i
       dphi=P_v02d(i,j-1)-P_v02d(i,j);
       %%%%oil
       if dphi >=0
           kro=Kro(Sw_v02d(i,j-1));
           uo=Uo(P_v02d(i,j-1));
           bo=Bo(P_v02d(i,j-1));
           Tf=kro/(uo*bo);
           Tx_o(i,j)=Tx_g(i,j)*Tf;
       else
           kro=Kro(Sw_v02d(i,j));
           uo=Uo(P_v02d(i,j));
           bo=Bo(P_v02d(i,j));
           Tf=kro/(uo*bo);
           Tx_o(i,j)=Tx_g(i,j)*Tf;
       end
       
       %%%water
        if dphi >=0
           krw=Krw(Sw_v02d(i,j-1));
           uw=1;
           bw=Bw(P_v02d(i,j-1));
           Tf=krw/(uw*bw);
           Tx_w(i,j)=Tx_g(i,j)*Tf;
       else
           krw=Krw(Sw_v02d(i,j));
           uw=1;
           bw=Bw(P_v02d(i,j));
           Tf=krw/(uw*bw);
           Tx_w(i,j)=Tx_g(i,j)*Tf;
       end
       
    end
end


%%%%% Y-direction

Ty_o=zeros(ny+1,nx);
Ty_w=zeros(ny+1,nx);
for i=2:ny
    for j=1:nx
        
        %%% potential of flow from j-1 to j
       dphi=P_v02d(i-1,j)-P_v02d(i,j);
       %%%%oil
       if dphi >=0
           kro=Kro(Sw_v02d(i-1,j));
           uo=Uo(P_v02d(i-1,j));
           bo=Bo(P_v02d(i-1,j));
           Tf=kro/(uo*bo);
           Ty_o(i,j)=Ty_g(i,j)*Tf;
       else
           kro=Kro(Sw_v02d(i,j));
           uo=Uo(P_v02d(i,j));
           bo=Bo(P_v02d(i,j));
           Tf=kro/(uo*bo);
           Ty_o(i,j)=Ty_g(i,j)*Tf;
       end
       
       %%%water
        if dphi >=0
           krw=Krw(Sw_v02d(i-1,j));
           uw=1;
           bw=Bw(P_v02d(i-1,j));
           Tf=krw/(uw*bw);
           Ty_w(i,j)=Ty_g(i,j)*Tf;
       else
           krw=Krw(Sw_v02d(i,j));
           uw=1;
           bw=Bw(P_v02d(i,j));
           Tf=krw/(uw*bw);
           Ty_w(i,j)=Ty_g(i,j)*Tf;
       end
       
    end
end
end


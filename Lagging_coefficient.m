
%%% grid and initial properties%%%%%%%%%%%%%%%%%%%%
grid.nx=15;
grid.ny=15;
grid.nz=1;
grid.Pinit=3000;
grid.Swinit=0.25;
grid.Soinit=0.75;
load('permx_final.mat');
grid.kx=reshape(permx,15,15);
grid.kx=grid.kx';
grid.ky=grid.kx;
grid.dx=30;
grid.dy=30;
grid.dz=30;

%%%%%%%%%property allocation

dx=grid.dx;
dy=grid.dy;
dz=grid.dz;
nx=grid.nx;
ny=grid.ny;
nz=grid.nz;
kx=grid.kx;
ky=grid.ky;
Pinit=grid.Pinit;
Swinit=grid.Swinit;
Soinit=grid.Soinit;
rw=0.35;
Pbhp=2900;
Qi=300*5.615;
s=+2.5;%%%%%%%%%skin
%%%%%%geometric transmissibilities

%%% X -direction
Tx_g=zeros(ny,nx+1);
for i=1:ny
    for j=2:nx
        K=2*kx(i,j-1)*kx(i,j)/(kx(i,j-1)+kx(i,j));
        Tx_g(i,j)=0.00633*K*dy*dz/dx;
    end
end

%%% Y-direction
Ty_g=zeros(ny+1,nx);
for i=2:ny
    for j=1:nx
        K=2*ky(i-1,j)*ky(i,j)/(ky(i-1,j)+ky(i,j));
        Ty_g(i,j)=0.00633*K*dx*dz/dy;
    end
end
       
%%%%%%initialize simulator

P_n=Pinit*ones(nx*ny,1);
Sw_n=Swinit*ones(nx*ny,1);
So_n=Soinit*ones(nx*ny,1);
PRESSURE(:,1)=P_n;
SATURATION(:,1)=Sw_n;
T_max=365;
dt=1;
t=0;
timestep=1;
 while t<T_max
    
    P_n2d=reshape(P_n,nx,ny);
    P_n2d=P_n2d';
    Sw_n2d=reshape(Sw_n,nx,ny);
    Sw_n2d=Sw_n2d';
    
    %%%%% calculation of fluid tansmissibilities%%%%%%%%%%%%
    
    %%% X -direction
Tx_o=zeros(ny,nx+1);
Tx_w=zeros(ny,nx+1);
for i=1:ny
    for j=2:nx
        
        %%% potential of flow from i-1 to i
       dphi=P_n2d(i,j-1)-P_n2d(i,j);
       %%%%oil
       if dphi >=0
           kro=Kro(Sw_n2d(i,j-1));
           uo=Uo(P_n2d(i,j-1));
           bo=Bo(P_n2d(i,j-1));
           Tf=kro/(uo*bo);
           Tx_o(i,j)=Tx_g(i,j)*Tf;
       else
           kro=Kro(Sw_n2d(i,j));
           uo=Uo(P_n2d(i,j));
           bo=Bo(P_n2d(i,j));
           Tf=kro/(uo*bo);
           Tx_o(i,j)=Tx_g(i,j)*Tf;
       end
       
       %%%water
        if dphi >=0
           krw=Krw(Sw_n2d(i,j-1));
           uw=1;
           bw=Bw(P_n2d(i,j-1));
           Tf=krw/(uw*bw);
           Tx_w(i,j)=Tx_g(i,j)*Tf;
       else
           krw=Krw(Sw_n2d(i,j));
           uw=1;
           bw=Bw(P_n2d(i,j));
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
       dphi=P_n2d(i-1,j)-P_n2d(i,j);
       %%%%oil
       if dphi >=0
           kro=Kro(Sw_n2d(i-1,j));
           uo=Uo(P_n2d(i-1,j));
           bo=Bo(P_n2d(i-1,j));
           Tf=kro/(uo*bo);
           Ty_o(i,j)=Ty_g(i,j)*Tf;
       else
           kro=Kro(Sw_n2d(i,j));
           uo=Uo(P_n2d(i,j));
           bo=Bo(P_n2d(i,j));
           Tf=kro/(uo*bo);
           Ty_o(i,j)=Ty_g(i,j)*Tf;
       end
       
       %%%water
        if dphi >=0
           krw=Krw(Sw_n2d(i-1,j));
           uw=1;
           bw=Bw(P_n2d(i-1,j));
           Tf=krw/(uw*bw);
           Ty_w(i,j)=Ty_g(i,j)*Tf;
       else
           krw=Krw(Sw_n2d(i,j));
           uw=1;
           bw=Bw(P_n2d(i,j));
           Tf=krw/(uw*bw);
           Ty_w(i,j)=Ty_g(i,j)*Tf;
       end
       
    end
end
        

%%%%% formulation of compressibility matrix
ctr=0;
D=zeros(2,2);
for i=1:ny
    for j=1:nx
        
        d11=0;
        d12=poro(P_n2d(i,j))*(1/Bw(P_n2d(i,j)));
        d21=0;
        d22=-poro(P_n2d(i,j))*(1/Bo(P_n2d(i,j)));
        
        ctr=ctr+1;
        d=[d11 d12;d21 d22];
        if ctr>1
         D=[D; d];
        else
         D=d;
        end
        
    end
end

   D= ((dx*dy*dz)/dt)*D;
   D=mat2cell(D,2*ones(nx*ny,1),2);
   D=blkdiag(D{:,1});
   
   %%%%%%%%%% formulation of transmissibility matrix
   %% center diagonal
   ctr=0;
   C=zeros(2,2);
   
   ro=0.2*dx;
   uw=1;
   for i=1:ny
       for j=1:nx
           a11=-(Ty_w(i+1,j)+Ty_w(i,j)+Tx_w(i,j+1)+Tx_w(i,j));
           a12=0;
           a21=-(Ty_o(i+1,j)+Ty_o(i,j)+Tx_o(i,j+1)+Tx_o(i,j));
           a22=0;
           
           if i==15 && j==15 %%% producer well
            
             WI=0.00633*(kx(i,j)*ky(i,j))^0.5*dz*(1/(log(ro/rw)+s));
             WIf_w=Krw(Sw_n2d(i,j))/(uw*Bw(P_n2d(i,j)));
             WIf_o=Kro(Sw_n2d(i,j))/(Uo(P_n2d(i,j))*Bo(P_n2d(i,j)));
             
             WI_w=WI*WIf_w;
             WI_o=WI*WIf_o;
             a11=a11-WI_w;
             a21=a21-WI_o;
           end
             
           a=[a11 a12;a21 a22];
           ctr=ctr+1;
           if ctr>1
               C=[C; a];
           else
               C=a;
           end
           
       end
   end
   
   C=mat2cell(C,2*ones(nx*ny,1),2);
   C=blkdiag(C{:,1});
   
   %%%%i+1 diagonal
   
   ctr=0;
   B=zeros(2,2);
   for i=1:ny
       for j=1:nx
           a11=Tx_w(i,j+1);
           a12=0;
           a21=Tx_o(i,j+1);
           a22=0;
           a=[a11 a12;a21 a22];
           ctr=ctr+1;
           if ctr>1
               B=[B ;a];
           else
               B=a;
           end
       end
   end
   B=mat2cell(B,2*ones(nx*ny,1),2);
   B{1,1}=[zeros(2,2),B{1,1}];
   B=blkdiag(B{1:end-1,1});
   B=[B;zeros(2,2*nx*ny)];
   
   %%%%i-1 diagonal
   
   ctr=0;
   A=zeros(2,2);
   for i=1:ny
       for j=1:nx
           a11=Tx_w(i,j);
           a12=0;
           a21=Tx_o(i,j);
           a22=0;
           a=[a11 a12;a21 a22];
           ctr=ctr+1;
           if ctr>1
               A=[A; a];
           else
               A=a;
           end
       end
   end
   A=mat2cell(A,2*ones(nx*ny,1),2);
   A=blkdiag(A{2:end,1});  
   A=[zeros(2,2*nx*ny-2);A];    
   A=[A,zeros(2*nx*ny,2)];
   
  %%%%j+1 diagonal
  
  ctr=0;
  E=zeros(2,2);
  for i=1:ny
      for j=1:nx
           a11=Ty_w(i+1,j);
           a12=0;
           a21=Ty_o(i+1,j);
           a22=0;
           a=[a11 a12;a21 a22];
           ctr=ctr+1;
           if ctr>1
               E=[E; a];
           else
               E=a;
           end
       end
  end  
     
  E=mat2cell(E,2*ones(nx*ny,1),2);
  E{1,1}=[zeros(2,2*nx),E{1,1}];
  E=blkdiag(E{1:end-nx,1});
  E=[E;zeros(2*nx,2*nx*ny)];
   %%%% j-1 diagonal
   
   ctr=0;
  F=zeros(2,2);
  for i=1:ny
      for j=1:nx
           a11=Ty_w(i,j);
           a12=0;
           a21=Ty_o(i,j);
           a22=0;
           a=[a11 a12;a21 a22];
           ctr=ctr+1;
           if ctr>1
               F=[F ;a];
           else
               F=a;
           end
       end
  end  
  
  F=mat2cell(F,2*ones(nx*ny,1),2);
  F=blkdiag(F{nx+1:end,1});
  F=[zeros(2*nx,2*nx*ny-2*nx);F];
  F=[F,zeros(2*nx*ny,2*nx)];
  
  T=E+B+C+A+F;
  
  
  %%%% formulation of rate matrix
  
 Q=zeros(2,1);
   ro=0.2*dx;
   uw=1;
   
  for i=1:ny
      for j=1:nx
          
         if i==15 && j==15%%%%% producer
            
             WI=0.00633*(kx(i,j)*ky(i,j))^0.5*dz*(1/(log(ro/rw)+s));
             WIf_w=Krw(Sw_n2d(i,j))/(uw*Bw(P_n2d(i,j)));
             WIf_o=Kro(Sw_n2d(i,j))/(Uo(P_n2d(i,j))*Bo(P_n2d(i,j)));
             
             WI_w=WI*WIf_w;
             WI_o=WI*WIf_o;
             
             qw=-WI_w*Pbhp;
             qo=-WI_o*Pbhp;
             
             q=[qw;qo];
             Q=[Q;q];
         elseif i==1 && j==1 %%%%% injector
             qw=-Qi;
             qo=0;
             q=[qw;qo];
             Q=q;
         else
             Q=[Q ;zeros(2,1)];
         end
      end
  end
  
  %%% formulating X_n matrix
  
 
  
  for i=1:ny
      for j=1:nx
          if i==1 && j==1
           X_n=[P_n2d(1,1);Sw_n2d(1,1)];
          else
          x=[P_n2d(i,j);Sw_n2d(i,j)];
          X_n=[X_n;x];
          end
      end
  end
  
  %%%%problem
  %% TX_n+1=D(X_n+1-X_n)+Q
  %%(T-D)X_n+1=-D*X_n+Q
  
 A_f=T-D;
 B_f=-D*X_n+Q;
 X_n1=A_f\B_f;
 
 %% update the values at n time level to be the new computed values
 
 ctr=0;
 
 for i=1:2:2*nx*ny-1
     ctr=ctr+1;
     P_n(ctr,1)=X_n1(i);
     Sw_n(ctr,1)=X_n1(i+1);
 end
 
 t=t+dt;%%%% increment time
 timestep=timestep+1;
 
 %%%% assigning the values to the matrix which stores all the pressures and
 PRESSURE(:,timestep)=P_n;
 SATURATION(:,timestep)=Sw_n;
 
 %%%% calculating rate
 
 i=15;%%% producer well grid
 j=15;
ro=0.2*dx; 
WI=0.00633*(kx(i,j)*ky(i,j))^0.5*dz*(1/(log(ro/rw)+s));
WIf_w=Krw(Sw_n(225,1))/(uw*Bw(P_n(225,1)));
WIf_o=Kro(Sw_n(225,1))/(Uo(P_n(225,1))*Bo(P_n(225,1)));
             
WI_w=WI*WIf_w;
WI_o=WI*WIf_o;

Rate_w(timestep-1)=WI_w*(P_n(225,1)-Pbhp)/5.615;
Rate_o(timestep-1)=WI_o*(P_n(225,1)-Pbhp)/5.615;
Time(timestep-1)=t;

end

%%%% cumulative production

for i=1:timestep-1
    if i==1
    Q_w(i)=Rate_w(i)*Time(i);
    Q_o(i)=Rate_o(i)*Time(i);
    else
    Q_w(i)=Q_w(i-1)+Rate_w(i)*Time(i);
    Q_o(i)=Q_o(i-1)+Rate_o(i)*Time(i);
    end
end

             
  
  
      
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FULLY IMPLICIT METHOD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
s=+2.5;%%%%skin
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

%%%%geometric well index
i=15;
j=15;
ro=0.2*dx;
WI_g=0.00633*(kx(i,j)*ky(i,j))^0.5*dz*(1/(log(ro/rw)+s));

%%%%%%initialize simulator

P_n=Pinit*ones(nx*ny,1);
Sw_n=Swinit*ones(nx*ny,1);
So_n=Soinit*ones(nx*ny,1);
T_max=365;

PRESSURE=zeros(nx*ny,T_max+1);
SATURATION=zeros(nx*ny,T_max+1);
PRESSURE(:,1)=P_n;
SATURATION(:,1)=Sw_n;

dt=1;
t=0;
timestep=1;
e_p=1e-06;
e_s=1e-03;
Niter=10;
Rate_w=zeros(1,T_max);
Rate_o=zeros(1,T_max);


%%%%%%%%%%%

while t<T_max
    
P_n1=1.05*P_n;
Sw_n1=1.005*Sw_n;

P_n2d=reshape(P_n,nx,ny);
P_n2d=P_n2d';
Sw_n2d=reshape(Sw_n,nx,ny);
Sw_n2d=Sw_n2d';

P_v0=P_n1;%%%%%%assign the pressure and saturation at zero level of iteration
Sw_v0=Sw_n1;
iter=0;
 disp('simulation time');
 disp(t);
while iter<Niter 
    %%% iteration%%
    P_v02d=reshape(P_v0,nx,ny);
    P_v02d=P_v02d';
    
    Sw_v02d=reshape(Sw_v0,nx,ny);
    Sw_v02d=Sw_v02d';
    
    %%%%%%Transmissibility calculation%%%%
    
    [ Tx_o,Ty_o,Tx_w,Ty_w ] = Trans(P_v02d,Sw_v02d,Tx_g,Ty_g);

%%%%%%%%%%residual vector calculation%%%%%
Res=zeros(2*nx*ny,2);
ctr=0;
for i=1:ny
    for j=1:nx
       Res_w = Rw( P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g);
       Res_o = Ro( P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g);
       ctr=ctr+1;
       Res(2*ctr-1,1)=Res_w;
       Res(2*ctr,1)=Res_o;
    end
end
Res=sparse(Res);

%%%%%%Jacobian calculation
e1=1e-06;
%%% centre diagonal%%%%
C=zeros(2*nx*ny,2);
ctr=0;
for i=1:ny
    for j=1:nx
        
        c11= Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'ij','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        c12= Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','ij',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        c21= Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'ij','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        c22= Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','ij',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        
        ctr=ctr+1;
        C(2*ctr-1,1)=c11/e1;
        C(2*ctr-1,2)=c12/e1;
        C(2*ctr,1)=c21/e1;
        C(2*ctr,2)=c22/e1;
 
    end
end

 C=mat2cell(C,2*ones(nx*ny,1),2);
 C=blkdiag(C{:,1});
 C=sparse(C);
   
B=zeros(2*nx*ny,2);
ctr=0;
%%%i+1 diagonal
for i=1:ny
    for j=1:nx
        
        if j==nx
        ctr=ctr+1;
        B(2*ctr-1,1)=0;
        B(2*ctr-1,2)=0;
        B(2*ctr,1)=0;
        B(2*ctr,2)=0;
        
        else
        b11=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'+i','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        b12=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','+i',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        b21=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'+i','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        b22=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','+i',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
       
        ctr=ctr+1;
        B(2*ctr-1,1)=b11/e1;
        B(2*ctr-1,2)=b12/e1;
        B(2*ctr,1)=b21/e1;
        B(2*ctr,2)=b22/e1;
        end
        
        
    end
end

   B=mat2cell(B,2*ones(nx*ny,1),2);
   B{1,1}=[zeros(2,2),B{1,1}];
   B=blkdiag(B{1:end-1,1});
   B=[B;zeros(2,2*nx*ny)];
   B=sparse(B);

  %%i-1 diagonal
  A=zeros(2*nx*ny,2);
  ctr=0;
  for i=1:ny
    for j=1:nx
        
        if j==1
         ctr=ctr+1;
        A(2*ctr-1,1)=0;
        A(2*ctr-1,2)=0;
        A(2*ctr,1)=0;
        A(2*ctr,2)=0;
        else
        a11=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'-i','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        a12=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','-i',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        a21=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'-i','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        a22=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','-i',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
       
        ctr=ctr+1;
        A(2*ctr-1,1)=a11/e1;
        A(2*ctr-1,2)=a12/e1;
        A(2*ctr,1)=a21/e1;
        A(2*ctr,2)=a22/e1;
        end
        
        
    end
  end
  
   A=mat2cell(A,2*ones(nx*ny,1),2);
   A=blkdiag(A{2:end,1});  
   A=[zeros(2,2*nx*ny-2);A];    
   A=[A,zeros(2*nx*ny,2)];
   A=sparse(A);
  %%%%%%j+1 diagonal
  
  E=zeros(2*nx*ny,2);
  ctr=0;
  for i=1:ny
    for j=1:nx
        
        if i==ny
        ctr=ctr+1;
        E(2*ctr-1,1)=0;
        E(2*ctr-1,2)=0;
        E(2*ctr,1)=0;
        E(2*ctr,2)=0;
        else
        e11=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'+j','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        e12=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','+j',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        e21=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'+j','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        e22=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','+j',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
       
        ctr=ctr+1;
        E(2*ctr-1,1)=e11/e1;
        E(2*ctr-1,2)=e12/e1;
        E(2*ctr,1)=e21/e1;
        E(2*ctr,2)=e22/e1;
        end
        
        
    end
  end

  E=mat2cell(E,2*ones(nx*ny,1),2);
  E{1,1}=[zeros(2,2*nx),E{1,1}];
  E=blkdiag(E{1:end-nx,1});
  E=[E;zeros(2*nx,2*nx*ny)];
  E=sparse(E);     
     %%%%%j-1 diagonal
     
     D=zeros(2*nx*ny,2);
     ctr=0;
    for i=1:ny
    for j=1:nx
        
        if i==1
        ctr=ctr+1;
        D(2*ctr-1,1)=0;
        D(2*ctr-1,2)=0;
        D(2*ctr,1)=0;
        D(2*ctr,2)=0;
        else
        d11=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'-j','none',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        d12=Rw(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','-j',Tx_w,Ty_w,Tx_g,Ty_g,dt,WI_g)-Res(2*j-1+2*(i-1)*nx);
        d21=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'-j','none',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
        d22=Ro(P_v02d,P_n2d,Sw_v02d,Sw_n2d,i,j,'none','-j',Tx_o,Ty_o,Tx_g,Ty_g,dt,WI_g)-Res(2*j+2*(i-1)*nx);
       
        ctr=ctr+1;
        D(2*ctr-1,1)=d11/e1;
        D(2*ctr-1,2)=d12/e1;
        D(2*ctr,1)=d21/e1;
        D(2*ctr,2)=d22/e1;
        end
        
    end
    end
  D=mat2cell(D,2*ones(nx*ny,1),2);
  D=blkdiag(D{nx+1:end,1});
  D=[zeros(2*nx,2*nx*ny-2*nx);D];
  D=[D,zeros(2*nx*ny,2*nx)];
  D=sparse(D);
  
  J= C+B+A+E+D;%%%%%JAcobian matrix
 J=sparse(J);
tic;
delta=-J\Res;
t2=toc;
% tic;
% tol=1e-10;
% setup.type = 'nofill'; % compute ILU(0)
% [L,U] = ilu(sparse(J), setup );
% maxit = 1000;
% [delta, flag, relres, iter_d] = bicgstab(sparse(J), Res, tol, maxit, L,U);
% t2=toc;
    
     ctr=0;
     deltaP_v1=zeros(nx*ny,1);
     deltaSw_v1=zeros(nx*ny,1);
 for i=1:2:2*nx*ny-1
     ctr=ctr+1;
     deltaP_v1(ctr,1)=delta(i,1);
     deltaSw_v1(ctr,1)=delta(i+1,1);
 end

     P_v1=zeros(nx*ny,1);
     Sw_v1=zeros(nx*ny,1);
 for i=1:225
     P_v1(i,1)=P_v0(i,1)+deltaP_v1(i,1);
     Sw_v1(i,1)=Sw_v0(i,1)+deltaSw_v1(i,1);
 end
 
 if norm(P_v1-P_v0)<=e_p && norm(Sw_v1-Sw_v0)<=e_s %%converge
    if max(Sw_v1)<=1 && min(Sw_v1)>=0%% saturation bounded
     P_n1=P_v1;
     Sw_n1=Sw_v1;
     break;
    end
 end
 
 if max(Sw_v1)>1 || min(Sw_v1)<0%%% saturation out of bound
     dt=dt/2;
     iter=0;
 end
 
 P_v0=P_v1;
 Sw_v0=Sw_v1;
 
 iter=iter+1;

 disp('iterations')
 disp(iter);
 
 if iter>Niter
     dt=dt/2;
     iter=0;
 end
 
 disp('timestep size');
 disp(dt);
end %%%%%%%%%%%%%%end of newton-raphson loop
 
 t=t+dt;
 timestep=timestep+1;
 PRESSURE(:,timestep)=P_n1;
 SATURATION(:,timestep)=Sw_n1;
 
 uw=1;
 WIf_w=Krw(Sw_n1(225,1))/(uw*Bw(P_n1(225,1)));
 WIf_o=Kro(Sw_n1(225,1))/(Uo(P_n1(225,1))*Bo(P_n1(225,1)));
 
 
 Rate_w(1,timestep-1)=WIf_w*WI_g*(P_n1(225,1)-Pbhp)/5.615;
 Rate_o(1,timestep-1)=WIf_o*WI_g*(P_n1(225,1)-Pbhp)/5.615;
 TIME(1,timestep-1)=dt;
 %%%Reassign pressure at n time level
 
 P_n=P_n1;
 Sw_n=Sw_n1;
     
 end %% end of timestep loop

 
    
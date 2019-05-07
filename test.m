nx=15;
ny=15;
i=0;
for t=30:30:360
    i=i+1;
    hfig=figure;
    image((reshape(SATURATION(:,t),nx,ny))','Cdatamapping','scaled');
    s=0;%%%skin
    S=num2str(s);
    T=num2str(t);
    Name=strcat('Water Saturation Map',',Skin-',S,',Time-',T,' days');
    title(Name);
    colorbar;
    F(i)=getframe(hfig);
    
end
movie(gcf,F,1,2);
movie2avi(F,'Saturation_LC_s0','fps',1);


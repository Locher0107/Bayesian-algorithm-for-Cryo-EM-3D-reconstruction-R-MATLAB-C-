Falcon_01=readMRCfile('Falcon_2012_06_12-17_17_05_0.mrc');


F = mat2gray(Falcon_01);% Use mat2gray to scale the image between 0 and 1
fig1=figure;
imshow(F);

saveas(fig1, '2Dmicrograph', 'jpg')
%%
%groel data set
a=readMRCfile('micrograph_001.mrcs');
a(:,:,1);
tt=mat2gray(a(:,:,2));
fig1=figure;
imshow(tt)
saveas(fig1, '2Dimage', 'jpg')
%%
Y=fft2(a(:,:,4));
shift_Y=fftshift(Y);
spectrum=log(1+abs(shift_Y));
plot(spectrum,'o')
imshow(mat2gray(abs(spectrum)))
%%
maximum=max(max(spectrum));
spectrum=255*spectrum/maximum;
spectrum=uint8(spectrum);
imshow(spectrum)



%%
[groel_coor, groel_type]=ReadPDBAtoms('4hel.pdb');
%%
plot3(groel_coor(1,:),  groel_coor(2,:), groel_coor(3,:), 'k'  );



%%
%%%%%%%%%%%%
%%%%%%%%%%%%%
%get CTF image
a=readMRCfile('micrograph_001.mrcs');
a(:,:,1);
tt=mat2gray(a(:,:,2));
imshow(tt)
% Angstrom per pixel= 1.42 A
%%
% to produce 2D CTF plot
lambda=0.020;
F=-23000;
Cs=2*10^7;
pix_a=1.42;

ctf2d=zeros(200, 200);
for( x=1:200)
    for (y=1:200)
        freq=sqrt((0-x)^2+(0-y)^2);
        freq=freq*1/(200*1.4);
        gamma_0=pi*lambda*F*freq^2+(1/2)*pi*Cs*lambda^3*freq^4;
        %2d ctf  value should be non-negative
        ctf2d(x,y)=abs( sin(gamma_0) );
        
        
    end
end
%%
fig1=figure;
imshow( fftshift(ctf2d));

%saveas(fig1, 'ctf_plot_3', 'jpg');



%%%
%new CTF 2D plot
%pixel size= 1.34
%%

%%
ctf2d=zeros(360, 360);
deltaf=23400;
Cs=2*10^7;
lambda=0.02;
K1=pi*lambda;
K2=(1/2)*pi*Cs*lambda^3;
Q0=0.1;
K3=sqrt(1-Q0^2);

defocus_average=-23400;
defocus_deviation=0;


xs=360*1.34;
ys=360*1.34;
for(i =0:359)
    ip=i;
    for(j=0:359)
        jp=j;
        
        x=jp/xs;
        y=ip/ys;
        
        u2=x^2+y^2;
        u=sqrt(u2);
        u4=u2*u2;

        ellipsoid_ang = atan2(y, x);
        cos_ellipsoid_ang_2 = cos(2*ellipsoid_ang);
        deltaf=defocus_average + defocus_deviation*cos_ellipsoid_ang_2;
        
        argument=K1*deltaf*u2+K2*u4;
        %K1 deltaf K2 K5
        
        if abs(argument)< (pi/2)
            retval=1;
        else
        retval=-(  K3*sin(argument)-Q0*cos(argument)  )   ;
        end
        
   
        ctf2d(i+1,j+1)= retval;
    end
end

imshow(ctf2d)
    

%%
%Test for spider file
ctf2d=zeros(360, 360);
NX=360;
NY=360;
LSM=NX+2;
SCX=2.0/NX;
SCY=2.0/NY;
ACR=atan(0.1/(1-0.1))
CS=2.0*10^7;
DZ=23400;
FMAXSPFREQ=0.373;

LAMBDA=0.02;
IYC=NY/2+1;
AZZ=0;

for K=1:NY
    KY=K-1;
    if KY>IYC
        KY=KY-NY;
    end
    for I=1:2:NX
        KX=(I-1)/2;
        AK=FMAXSPFREQ* sqrt( (KX*SCX)^2+(KY*SCY)^2   );
        if (KX==0)
            AZ=pi/2;
        else
            AZ=atan(KY/KX)+pi/2;
        end
        AZR=AZZ*(pi/180);
        
        DZA=0;
        DZZ=DZ+DZA/2*sin(2*(AZ-AZR));
        
        %then call TFD
        F1=1/sqrt(CS*LAMBDA);
        F2=sqrt( sqrt(CS*LAMBDA^3) );
        %Q1=(Q*F2)^2;
        AKK=AK*F2; %AKK is the frequency
                            
        
        DZ1=DZZ*F1;%DZZ is the defocus; F1, F2 are some combinations of Cs and lambda
        
        QQT=2*pi*(0.25*AKK^4-0.5*DZ1*AKK^2);
        WGH=ACR;
        B=-sin(QQT-WGH);

        ctf2d( KX+1, K)=B;
        ctf2d(NX-KX, K)=B;
    end
end

%%
imshow(ctf2d)
%%
fig1=figure;
imshow(fftshift(ctf2d))

saveas(fig1, 'ctf_plot_3', 'jpg');
    

    







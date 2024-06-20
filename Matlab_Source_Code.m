%The Developed Matlab Source Code
t1=[0:0.01e-9:200e-9];
inter=0.01e-9; Fs=1/inter;
N=2^nextpow2(length(t1));
t2=(0:(N-1))*inter;
V1=0; V2=1;
TC1=3e-9; TD1=5e-9;
TC2=5e-9; TD2=25e-9;
voltage=zeros(1,length(t2));
for i=1:1:length(t2)
    if t2(i)<TD1;
        voltage(i)=V1;
    elseif t2(i)>=TD1&&t2(i)<=TD2;
        voltage(i)=V1+(V2-V1).*(1-exp(-(t2(i)-TD1)/TC1));
    else t2(i)>TD2;
        voltage(i)=V1+(V2-V1).*(1-exp(-(t2(i)-TD1)/TC1))+(V1-V2)*(1-exp(-(t2(i)-TD2)/TC2));
    end
end
plot(t2,voltage,'b','LineWidth',1.5)


Vs=fft(voltage,N);
f=Fs/N*(0:1:N/2);
r1=0.01; l1=15e-9; c1=15e-12; g1=0.001;
r2=0.01*5; l2=5e-9*5; c2=75e-12*5; g2=0.001*5;
r3=0.01*10; l3=15e-9*10; c3=25e-12*10; g3=0.001*10;
r0=0.02/2; l0=25e-9/2; c0=50e-12/2; g0=0.001/2;
Zs=0.1; L1=1; L2=2;
VAL=zeros(1,N/2+1); VBL=zeros(1,N/2+1); VCL=zeros(1,N/2+1);
for k=1:1:N/2+1
    z1=r1+1i.*2.*pi.*f(1,k).*l1; y1=g1+1i.*2.*pi.*f(1,k).*c1; Zc1=sqrt(z1./y1); Yc1=1./Zc1; gamma1=sqrt(z1.*y1);
    z2=r2+1i.*2.*pi.*f(1,k).*l2; y2=g2+1i.*2.*pi.*f(1,k).*c2; Zc2=sqrt(z2./y2); Yc2=1./Zc2; gamma2=sqrt(z2.*y2);
    z3=r3+1i.*2.*pi.*f(1,k).*l3; y3=g3+1i.*2.*pi.*f(1,k).*c3; Zc3=sqrt(z3./y3); Yc3=1./Zc3; gamma3=sqrt(z3.*y3);
    Ys1=Yc1./tanh(gamma1.*L2); Ym1=-Yc1./sinh(gamma1.*L2); Ys2=Yc2./tanh(gamma2.*L2); Ym2=-Yc2./sinh(gamma2.*L2); Ys3=Yc3./tanh(gamma3.*L2); Ym3=-Yc3./sinh(gamma3.*L2);
    W1=1; D1=Ys1; E1=-Ym1; F1=Ys1; G1=-Ym1; W2=1; D2=Ys2; E2=-Ym2; F2=Ys2; G2=-Ym2; W3=1; D3=Ys3; E3=-Ym3; F3=Ys3; G3=-Ym3;
    z0=r0+1i.*2.*pi.*f(1,k).*l0; y0=g0+1i.*2.*pi.*f(1,k).*c0; Zc0=sqrt(z0./y0); Yc0=1./Zc0; gamma0=sqrt(z0.*y0);
    K1=cosh(gamma0.*L1)+D1./W1.*Zc0.*sinh(gamma0.*L1)+F3./W3.*Zc0.*sinh(gamma0.*L1)+G1./W1.*Zc0.*sinh(gamma0.*L1);
    K2=-cosh(gamma0.*L1)-E1./W1.*Zc0.*sinh(gamma0.*L1)-D2./W2.*Zc0.*sinh(gamma0.*L1)-F1./W1.*Zc0.*sinh(gamma0.*L1);
    K3=-G3./W3.*Zc0.*sinh(gamma0.*L1)+E2./W2.*Zc0.*sinh(gamma0.*L1);
    K4=cosh(gamma0.*L1)+D1./W1.*Zc0.*sinh(gamma0.*L1)+F3./W3.*Zc0.*sinh(gamma0.*L1)+E3./W3.*Zc0.*sinh(gamma0.*L1);
    K5=-E1./W1.*Zc0.*sinh(gamma0.*L1)+G2./W2.*Zc0.*sinh(gamma0.*L1);
    K6=-cosh(gamma0.*L1)-G3./W3.*Zc0.*sinh(gamma0.*L1)-D3./W3.*Zc0.*sinh(gamma0.*L1)-F2./W2.*Zc0.*sinh(gamma0.*L1);
    K7=(K1.*K6-K3.*K4)./(K3.*K5-K2.*K6); K8=(K2.*K4-K1.*K5)./(K3.*K5-K2.*K6);
    B1=Zs.*Yc0.*sinh(gamma0.*L1)+cosh(gamma0.*L1)./3; B2=Zs.*cosh(gamma0.*L1)+Zc0.*sinh(gamma0.*L1)./3;
    Y1=D1./W1+F3./W3-G1./W1-E3./W3; Y2=D2./W2+F1./W1-E1./W1-G2./W2; Y3=D3./W3+F2./W2-G3./W3-E2./W2;
    VAL(1,k)=Vs(1,k)./(B1+Y1.*B2+K7.*(B1+Y2.*B2)+K8.*(B1+Y3.*B2));
    VBL(1,k)=K7.*VAL(1,k);
    VCL(1,k)=K8.*VAL(1,k);
end
VAL_1(1,1:1:N/2+1)=VAL(1,1:1:N/2+1);      VBL_1(1,1:1:N/2+1)=VBL(1,1:1:N/2+1);      VCL_1(1,1:1:N/2+1)=VCL(1,1:1:N/2+1);
VAL_1(1,N/2+2:1:N)=conj(VAL(1,N/2:-1:2)); VBL_1(1,N/2+2:1:N)=conj(VBL(1,N/2:-1:2)); VCL_1(1,N/2+2:1:N)=conj(VCL(1,N/2:-1:2));
VAL_2=ifft(VAL_1);                        VBL_2=ifft(VBL_1);                        VCL_2=ifft(VCL_1);
VAL_t=real(VAL_2);                        VBL_t=real(VBL_2);                        VCL_t=real(VCL_2);

subplot(1,3,1)
plot(t2,VAL_t,'r','LineWidth',1.5)
% hold on
% plot(pspice(:,1),pspice(:,2),'k--','LineWidth',1.5)
grid minor
xlim([0e-7 2e-7])

subplot(1,3,2)
plot(t2,VBL_t,'b','LineWidth',1.5)
% hold on
% plot(pspice(:,1),pspice(:,3),'m--','LineWidth',1.5)
grid minor
xlim([0e-7 2e-7])

subplot(1,3,3)
plot(t2,VCL_t,'k','LineWidth',1.5)
% hold on
% plot(pspice(:,1),pspice(:,4),'c--','LineWidth',1.5)
grid minor
xlim([0e-7 2e-7])
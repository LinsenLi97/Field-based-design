%The code for Fig. S4 in the supplemetary information of https://arxiv.org/abs/2101.02366
% Linsen Li
n0=17; %n0 slot for the simulation
W=[42 ,25  , 70 , 73  , 96 , 80  , 82  ,54  , 50  , 58 , 64  ,88  ,83  , 89 , 115 , 120 , 120];
xh=[0.2920 ,  0.5260  ,  0.7220  , 0.9200  , 1.1220  ,1.3400 , 1.5360 , 1.8880 , 2.2000 , 2.5180 , 2.8420 ,  3.1660  ,3.4960  , 3.8120 , 4.1560  ,  4.5140  ,  4.8940];
deltax=zeros(1,n0); % get the distance between different slots
deltax(1)=round((xh(1)-W(1)*1e-3/2)/0.002); % the input parameters will be integer.
for i=2:n0
    deltax(i)=round((xh(i)-W(i)*1e-3/2-W(i-1)*1e-3/2-xh(i-1))/0.002);
end
xxx=[W,deltax]; % get the general in put of for the transfer matrix model optimization. It will be 2n0 integer.

r0=zeros(1,10);
k0=20.12; %wave vector for 637nm slab mode in 150nm diamond membrane with refractive index n=2.41
deltaphi=k0*2e-3; %unit phase difference
m0=round(3.1415926*500/k0); % minimum pi phase distance
dmin=110; % minimum slot to slot distance
a0=0.23; % The target near field width (got from target far field)
W=round(xxx(1:n0));
W(W>120)=120; % max slot width 120nm
W(W<10)=10; % min slot width 10nm
deltax=round(xxx(n0+1:2*n0));
while (deltax(1)>298) %slot distance max
    deltax(1)=deltax(1)-m0;
end
while (deltax(1)<dmin) %slot distance min
    deltax(1)=deltax(1)+m0;
end
for (i=2:n0) %slot distance max for i = 2 ... n0 slot
while (deltax(i)>210)
    deltax(i)=deltax(i)-m0;
end
end
% look up table for reflection and transmission
rc=load("rc.mat");
rc=rc.rc;
tc=load("tc.mat");
tc=tc.tc;
%look up table for scattering electric field and magnetic field
abs0=load("abs0.mat");
abs0=abs0.abs0;
psi0=load("psi0.mat");
psi0=psi0.psi0;
absH=load("Hx_abs.mat");
absH=absH.Hx_abs';
psiH=load("Hx_angle.mat");
psiH=psiH.Hx_angle';
% transfer matrix calculation start
M0=zeros(2,2,n0);
W0=round(W/2); %index for searching the look up table (0:2:120)nm
for i=1:n0
    j=n0+1-i;
    M0(:,:,j)=[1/(tc(W0(j))) -rc(W0(j))/tc(W0(j)); (rc(W0(j))/tc(W0(j))) (tc(W0(j))^2-rc(W0(j))^2)/tc(W0(j))];
end
XR=zeros(2,n0);
XL=zeros(2,n0);
XR(1,n0)=1;
XR(2,n0)=0;
for i=1:n0-1
    XL(:,n0+1-i)=M0(:,:,n0+1-i)*XR(:,n0+1-i);
    XR(:,n0-i)=[exp(-1i*deltaphi*deltax(n0+1-i)),0;0,exp(1i*deltaphi*deltax(n0+1-i))]*XL(:,n0+1-i);
end
XL(:,1)=M0(:,:,1)*XR(:,1);
%Phase condition for the 1st slot location (left and right of the origin can have constructive inteference)
deltax(1)=round(((angle(XL(1,1))-angle(XL(2,1)))/(2*deltaphi)))+round(3.1415926/deltaphi); 
while (deltax(1)>299)
    deltax(1)=deltax(1)-m0;
end
while (deltax(1)<dmin)
    deltax(1)=deltax(1)+m0;
end
% Coherent adding the contribution
% The reason having four terms can consider the following figure: https://drive.google.com/file/d/1ndQFEYM4suqzr3VUH6Ej9XkUVoDuY2C7/view?usp=sharing
Ey=zeros(1,2501);
Hx=zeros(1,2501);
Eyf=zeros(1,5001);
Hxf=zeros(1,5001);
for i=1:n0
    a1=fliplr(abs0(W0(i),2501-xslot(i):5001-xslot(i)));
    p1=fliplr(psi0(W0(i),2501-xslot(i):5001-xslot(i)));   
    a2=abs0(W0(i),5001-xslot(i):7501-xslot(i));
    p2=psi0(W0(i),5001-xslot(i):7501-xslot(i));   
    a3=fliplr(abs0(W0(i),2501+xslot(i):5001+xslot(i)));
    p3=fliplr(psi0(W0(i),2501+xslot(i):5001+xslot(i)));   
    a4=abs0(W0(i),5001+xslot(i):7501+xslot(i));
    p4=psi0(W0(i),5001+xslot(i):7501+xslot(i));
    Ey=Ey+abs(XL(1,i))*a1.*exp(1i*(p1+ones(1,2501).*angle(XL(1,i))));
    Ey=Ey+abs(XL(1,i))*a2.*exp(1i*(p2+ones(1,2501).*angle(XL(1,i))));
    Ey=Ey+abs(XR(2,i))*a3.*exp(1i*(p3+ones(1,2501).*angle(XR(2,i))));
    Ey=Ey+abs(XR(2,i))*a4.*exp(1i*(p4+ones(1,2501).*angle(XR(2,i))));
% H field coherent adding 
    aH1=fliplr(absH(W0(i),2501-xslot(i):5001-xslot(i)));
    pH1=fliplr(psiH(W0(i),2501-xslot(i):5001-xslot(i)));   
    aH2=absH(W0(i),5001-xslot(i):7501-xslot(i));
    pH2=psiH(W0(i),5001-xslot(i):7501-xslot(i));   
    aH3=fliplr(absH(W0(i),2501+xslot(i):5001+xslot(i)));
    pH3=fliplr(psiH(W0(i),2501+xslot(i):5001+xslot(i)));   
    aH4=absH(W0(i),5001+xslot(i):7501+xslot(i));
    pH4=psiH(W0(i),5001+xslot(i):7501+xslot(i));
    Hx=Hx+abs(XL(1,i))*aH1.*exp(1i*(pH1+ones(1,2501).*angle(XL(1,i))));
    Hx=Hx+abs(XL(1,i))*aH2.*exp(1i*(pH2+ones(1,2501).*angle(XL(1,i))));
    Hx=Hx+abs(XR(2,i))*aH3.*exp(1i*(pH3+ones(1,2501).*angle(XR(2,i))));
    Hx=Hx+abs(XR(2,i))*aH4.*exp(1i*(pH4+ones(1,2501).*angle(XR(2,i))));
end
Ey_tar1D=zeros(1,5000);
for j=1:5001
    Ey_tar1D(j) = (exp(-(j-2501).^2/(0.8706*500)^2));
end
% result plot
figure(1)
plot(0:0.002:5,abs(Ey)/max(abs(Ey)),'r','linewidth',2);
hold on
plot(0:0.002:5,Ey_tar1D(2501:5001)/max(Ey_tar1D),'c','linewidth',2)
Tz=(conj(Ey(1:2001)).*(Hx(1:2001)));
Pz=-sum(real(Tz(1:2001)))*2*0.002;
Px=0.7126*1e-3;
P_ratio=Pz/(2*(Px+Pz));
xlabel('x(um)')
ylabel('Ey')
legend("TMM","Target")
modeoverlap3=sum((Ey(1:2501)).*Ey_tar1D(2501:5001))/sqrt(sum(abs(Ey(1:2501)).^2)*sum(abs(Ey_tar1D(2501:5001)).^2));
modeoverlap=abs(modeoverlap3)
FOM=modeoverlap4*modeoverlap4*P_ratio


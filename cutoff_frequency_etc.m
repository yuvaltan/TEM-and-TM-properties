close all
clear all
clc

global epsilon1 epsilon2 epsilon3 miu dist z1 z2 z3 N f0 w0 d tau
d=2*10^-3;
miu=4*pi*10^-7;
dist=2*10^-3;
N=1024
epsilon1=3.9*8.85*10^-12;
epsilon2=sqrt(3.9*2.5)*8.85*10^-12;
epsilon3=2.5*8.85*10^-12;
f0= 40*10^9;
w0=2*pi*f0;
beta2=2*pi*f0*sqrt(miu*epsilon2);
lamda2=2*pi/beta2;
z1=sqrt(miu/epsilon1)*2/10;
z2=sqrt(miu/epsilon2)*2/10;
z3=sqrt(miu/epsilon3)*2/10;
l2=lamda2/4;
v2=1/sqrt(miu*epsilon2);
w=linspace(-4*w0,4*w0,1024);
tau = 60*10^-12;

%---------------------
% smith plot
%---------------------
ret_co_l2=(z3-z2)/(z3+z2) ;
beta2_w= sqrt(miu*epsilon2)*w;
z_in = z2*(1+ ret_co_l2*exp(-2*1i*l2*beta2_w))./(1-ret_co_l2*exp(-2*1i*l2*beta2_w));
gamma_in = (z_in-z1)./(z_in+z1);
smithplot(gamma_in);

%----------------------
% transfer coeffcient
% ---------------------
ret_co_l1=(((1+ret_co_l2*exp(-2i*beta2*l2))*z2)-((1-ret_co_l2*exp(-2i*beta2*l2))*z1))/(((1+ret_co_l2*exp(-2i*beta2*l2))*z2)+((1-ret_co_l2*exp(-2i*beta2*l2))*z1));
trans_co3=@(w) (1+ret_co_l1)*(1+ret_co_l2)*(exp(-1i.*w.*sqrt(miu*epsilon2)*l2))./(1+(ret_co_l2*(exp(-2i.*w.*sqrt(miu*epsilon2)*l2))));
v0 =@(w) (0.5*sqrt(pi)*tau*(exp(-tau^2/4.*(w-w0).^2)+exp(-tau^2/4.*(w+w0).^2)));%result ç
v3p = @(w) trans_co3(w).* v0(w); %result 

figure
plot(w/w0,abs(trans_co3(w)))
title('transfer coefficient line 3|');
xlabel('w/w0');
ylabel('|tau_3|');
grid on

%---------------------------------------------
% bandwidth for v0 and v3
%---------------------------------------------
i=0;
for q=-4*w0:w0/128:4*w0
    i=i+1;
   w_v(i)=q;% omega indexing
   v3_v(i)=v3p(q);%sampling v3 in position z=0
end
M_v3=max(abs(v3_v));
[a,b]=max(abs(v3_v));
pos_w=[0,0];
stoper=0;
for j = 1:length(v3_v)
     
    if abs(v3_v(j))>=M_v3*exp(-1) && stoper==0
        pos_w(1)=j;%start of bandwith
        stoper=1;
    elseif abs(v3_v(j))<=M_v3*exp(-1) && stoper==1
        pos_w(2)=j;%end of bandwith
        stoper=2;
        break;
    end
end
bandwidth = abs(w_v(pos_w(1))-w_v(pos_w(2)))/(2*pi);

% ---------------------------------------------
% abs|v0(w)|
% ---------------------------------------------
figure
plot(w/w0,abs(v0(w)))
title('|v0(w)| normalized graph of amplitude as functin of w');
xlabel('w/w0');
ylabel('|v0(w)|');
grid on

figure
plot(w/w0,abs(v3p(w)))
title('|v3(w)|normalized graph of amplitude as functin of w');
xlabel('w/w0');
ylabel('|v3(w)|');
grid on

%---------------------------------------------
% inverse Fourie transform of v3
%---------------------------------------------
dt=pi/(4*w0);
samp_v3=circshift(v3p(w),N/2);% samples of v3
fourie_v3= ifft(samp_v3)/dt;
inverse_fourier_v3=circshift(fourie_v3,N/2);%back to tim÷
time_div = -1023*dt/2:dt:1023*dt/2;
g_t = exp(-(time_div-l2/v2).*(time_div-l2/v2)/tau^2);%the function g from the begining

figure
plot(time_div*(10^9),inverse_fourier_v3)
hold on
plot(time_div*(10^9), g_t ,'g')
title('v3(t)');
xlabel('t(10^-9)');
ylabel('v3(t)');
grid on
legend('v3(t)','g(t)');
g_t = exp(-(time_div-l2/v2).*(time_div-l2/v2)/tau^2);    %the function g from the begining


%----------------------
% beta_n
%----------------------
wc=pi/(d*sqrt(miu*epsilon1));
beta=@(n,w) (double((1*w>wc)-1*(w<wc))).*sqrt(w.^2*miu*epsilon1-(n.*pi/d)^2); %beta_n(w)

figure   %for the first three betas
plot(w,real(beta(1,w)))
hold on
plot(w,real(beta(2,w)))
plot(w,real(beta(3,w)))
title('Real[beta(w)]');
xlabel('w/w0');
ylabel('Re{beta(w)}');
grid on
legend('beta1(w)','beta2(w)','beta3(w)');

%--------------------------
% v1(w) for different z 
%--------------------------

v1=@(w,z) v0(w).*exp(-1i*z.*beta(1,w));   % v for basic character (n=1)
 
figure
plot(w/w0,log(abs(v1(w,0))))  %z=0
hold on
plot(w/w0,log(abs(v1(w,5*10^-3))))  %z=0.5*10^2
plot(w/w0,log(abs(v1(w,10*10^-3))))  %z=10^2
title('log|v1(w)| for different z[mm]');
xlabel('w/w0');
ylabel('log|v1(w)|');
grid on
legend('v1(z=0,w)','v1(z=5,w)','v1(z=10,w)');
 
figure
plot(w/w0,angle(v1(w,0)))
hold on
plot(w/w0,angle(v1(w,5*10^-3)))
plot(w/w0,angle(v1(w,10*10^-3)))
title('phase of v1(w) for different z[mm])');
xlabel('w/w0');
ylabel('phase(v1(w))');
grid on
legend('v1(z=0,w)','v1(z=5,w)','v1(z=10,w)');
 

%---------------------------------------------
% inverse fourier for each z
%---------------------------------------------
samp_v1_0=circshift(v1(w,0),512);    %fourier for each z
fourier_v1_0= ifft(samp_v1_0)/dt;
inverse_v1_0=circshift(fourier_v1_0,512);

samp_v1_5=circshift(v1(w,5*10^-3),512);
fourier_v1_5= ifft(samp_v1_5)/dt;
inverse_v1_5=circshift(fourier_v1_5,512);

samp_v1_10=circshift(v1(w,10*10^-3),512);
fourier_v1_10= ifft(samp_v1_10)/dt;
inverse_v1_10=circshift(fourier_v1_10,512);

figure
plot(time_div*10^9,inverse_v1_0)
hold on
plot(time_div*10^9,inverse_v1_5)
plot(time_div*10^9,inverse_v1_10)
plot(time_div*10^9,g_t)
title('v1(t) for different z');
xlabel('t[10^-9]');
ylabel('v1(t)');
grid on
legend('v1(z=0,t)','v1(z=5,t)','v1(z=10,t)','g(t)');
 

%---------------------------------------------
% taylor of beta for each z
%---------------------------------------------
beta_0 = beta(1,w0);    %zero order taylor coe
beta_1 = (w0*miu*epsilon1)/sqrt((w0^2*miu*epsilon1-(pi/d)^2));    %first order talor coe
v1new_taylor=@(z) exp(-(time_div-beta_1*z).*(time_div-beta_1*z)/tau^2).*cos(w0*time_div-beta_0*z); 

figure
plot(time_div*10^9,inverse_v1_0)
hold on
plot(time_div*10^9,v1new_taylor(0))
title('v1(t) for z=0mm');
xlabel('t(10^-9)');
ylabel('v1(t)');
grid on
legend('v1(t)','v1(t)-taylor');
 
figure
plot(time_div*10^9,inverse_v1_5)
hold on
plot(time_div*10^9,v1new_taylor(5*10^-3))
title('v1(t) for z=5mm');
xlabel('t(10^-9)');
ylabel('v1(t)');
grid on
legend('v1(t)','v1(t)-taylor');

figure
plot(time_div*10^9,inverse_v1_10)
hold on
plot(time_div*10^9,v1new_taylor(10^-2))
title('v1(t) for z=10mm');
xlabel('t(10^-9)');
ylabel('v1(t)');
grid on
legend('v1(t)','v1(t)-taylor');
 

% ---------------------------------------------
% q13 q14 for tau=600ps 
% ---------------------------------------------
tau=600*10^-12;  
v0_new =@(w) (0.5*sqrt(pi)*tau*(exp(-tau^2/4.*(w-w0).^2)+exp(-tau^2/4.*(w+w0).^2)));
v1_new=@(w,z) v0_new(w).*exp(-1i*z.*beta(1,w));
g_t = exp(-(time_div-l2/v2).*(time_div-l2/v2)/tau^2);

samp_v1_0=circshift(v1_new(w,0),512);    %fourier for each z
fourier_v1_0= ifft(samp_v1_0)/dt;
inverse_v1_0=circshift(fourier_v1_0,512);

samp_v1_5=circshift(v1_new(w,5*10^-3),512);
fourier_v1_5= ifft(samp_v1_5)/dt;
inverse_v1_5=circshift(fourier_v1_5,512);

samp_v1_10=circshift(v1_new(w,10*10^-3),512);
fourier_v1_10= ifft(samp_v1_10)/dt;
inverse_v1_10=circshift(fourier_v1_10,512);

figure
plot(time_div*10^9,inverse_v1_0)
hold on
plot(time_div*10^9,inverse_v1_5)
plot(time_div*10^9,inverse_v1_10)
plot(time_div*10^9,g_t)
title('v1(t) for different z');
xlabel('t[10^-9]');
ylabel('v1(t)');
grid on
legend('v1(z=0,t)','v1(z=5,t)','v1(z=10,t)','g(t)');
 



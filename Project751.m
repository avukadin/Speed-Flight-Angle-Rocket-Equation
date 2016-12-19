clc;
clear variables;
beep off;
figure('Color',[1 1 1]);

%====================================
%Paramters
%====================================
h=.1;%Time step size
T=2000;%Final time
v=zeros(round(T/h)+1,1);%Speed
theta=zeros(round(T/h)+1,1);%Flight angle
v(1)=4;%Initial speed
theta(1)=pi/2.2;%Initial angle
u= 1000;%Magnitude of thrust vector 
Fg=@(t) 9.8; %Gravitational acceleration
tau=10; %Burn rate of fuel
MI =500; %Initial mass of rocket
MF=MI*(0.8); %Initial mass of fuel
ME=MI-MF; %Mass of empty rocket
M=@(t) ME+MF*exp(-t/tau); %Mass of rocket at time t
Mp=@(t) -MF*exp(-t/tau)/tau; %Derivative of M(t)
c=@(t) 10; %Drag

f=@(t,v,theta)-Mp(t)*u/M(t) -Fg(t).*sin(theta)-c(t)*(v)/M(t);
g=@(t,v,theta) -Fg(t).*cos(theta)/(v);

%====================================
%4th Order Runge-Kutta Discretization
%====================================
for i=2:size(v,1)
    
    k0=h*f(h*(i-1),v(i-1),theta(i-1));
    l0=h*g(h*(i-1),v(i-1),theta(i-1));
    
    k1=h*f(h*(i-1)+0.5*h,v(i-1)+0.5*k0,theta(i-1)+0.5*l0);
    l1=h*g(h*(i-1)+0.5*h,v(i-1)+0.5*k0,theta(i-1)+0.5*l0);
    
    k2=h*f(h*(i-1)+0.5*h,v(i-1)+0.5*k1,theta(i-1)+0.5*l1);
    l2=h*g(h*(i-1)+0.5*h,v(i-1)+0.5*k1,theta(i-1)+0.5*l1);
    
    k3=h*f(h*(i-1)+h,v(i-1)+k2,theta(i-1)+l2);
    l3=h*g(h*(i-1)+h,v(i-1)+k2,theta(i-1)+l2);
    
    v(i)=v(i-1)+(1/6)*(k0+2*k1+2*k2+k3);
    theta(i)=theta(i-1)+(1/6)*(l0+2*l1+2*l2+l3);
    
end
%====================================
%Plotting of Orbits
%====================================
hold on
plot(v,theta,'LineWidth',2);
ylabel('Theta (Radians)') % x-axis label
xlabel('Speed (m/s)') % y-axis label




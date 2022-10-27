function dxdt = integratingfunc_Qn1_JeffcottRotor(t,x)
%INTEGRATINGFUNCTION Summary of this function goes here
%   Detailed explanation goes here

%Constants
g = 9.800; %m/s2
span = 1;    %m
dia = 0.025;   %m
m = 5;  %kg
E = 210e9;  %N/m2
I = (pi/64)*(dia^4);
k = (48*E*I)/(span^3);
k_e = 0.95*k;
k_t = 0.95*k;
w_not = (k/m)^0.5;
w_not_e = (k_e/m)^0.5;
w_not_t = (k_t/m)^0.5;
del = m*g/k;    
del_e = m*g/k_e;
del_t = m*g/k_t;
c = 100; %Ns/m
beta = pi/4;
eccen = 0.005; %m

%Dimensionless terms
theta_accn = 20; %rad/sec2
lambda = theta_accn/(w_not^2);
lambda_e = theta_accn/(w_not_e^2);
lambda_t = theta_accn/(w_not_t^2);
tau = w_not*t;
tau_e = w_not_e*t;
tau_t = w_not_t*t;
omega = lambda*tau;
omega_e = lambda_e*tau_e;
omega_t = lambda_t*tau_t;
theta = 0.5*lambda*(tau^2);
theta_e = 0.5*lambda_e*(tau_e^2);
theta_t = 0.5*lambda_t*(tau_t^2);
eta = c/(2*m*w_not);
eta_e = c/(2*m*w_not_e);
eta_t = c/(2*m*w_not_t);
e = eccen/del;
e_e = eccen/del_e;
e_t = eccen/del_t;

%Forces
Fe = (e_e*(omega_e^2)*cos(beta))+(e_e*lambda_e*sin(beta))-cos(theta_e);
Feta = (e_t*(omega_t^2)*sin(beta))-(e_t*lambda_t*cos(beta))+sin(theta_t);

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

dxdt = zeros(4,1);
dxdt(1) = x2;
dxdt(2) = Fe + (x3*lambda_e) + (2*x4*omega_e) + (x1*(omega_e^2)) - (2*eta_e*x2) + (2*eta_e*omega_e*x3) - x1;
dxdt(3) = x4;
dxdt(4) = Feta - (x1*lambda_t) - (2*x2*omega_t) + (x3*(omega_t^2)) - (2*eta_t*x4) - (2*eta_t*omega_t*x1) - x3;

end

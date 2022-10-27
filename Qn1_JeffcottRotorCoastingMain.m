tstart = 0;
tend = 50;
n = 100000;

%Constants

g = 9.800; %m/s2
span = 1;    %m
dia = 0.025;   %m
m = 5;  %kg
E = 210e9;  %N/m2
I = (pi/64)*(dia^4);
k = (48*E*I)/(span^3);
k_e = 0.9*k;
k_t = 0.97*k;
w_not = (k/m)^0.5;
w_not_e = (k_e/m)^0.5;
w_not_t = (k_t/m)^0.5;
del = m*g/k;    
del_e = m*g/k_e;
del_t = m*g/k_t;
c = 100; %Ns/m
beta = pi/4;
eccen = 0.005; %m

tspan = linspace(tstart,tend,n);
X1 = 0;
X2 = 0;
X3 = 0;
X4 = 0;
xinit = [X1,X2,X3,X4];
%options = odeset('RelTol',1e-3,'AbsTol',1e-3);
[t,x] = ode45(@integratingfunc_Qn1_JeffcottRotor, tspan, xinit);

plot(tspan,x(:,1)*del_e,tspan,x(:,3)*del_t);
xlabel('time(sec)');
ylabel('disp(m)');
legend('vertical','horizontal');
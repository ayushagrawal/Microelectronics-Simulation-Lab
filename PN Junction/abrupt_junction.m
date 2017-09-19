%% THE TECHNIQUES INVOLVED
%  space_charge_density, rho(x) = e.(h(x) - n(x) + N_D(x) - N_A(x))
%  Poisson Equation, epsilon . (d^2V(x)/dx^2) = -rho(x)

%% IMPLEMENTING THE TECHNIQUES

% constants
epsilon0 = 8.854 * 10^-12;
epsilon_si = 11.68;
epsilon = epsilon0*epsilon_si;
nm = 1e-9;      % nano meter -> m


k = 1.38e-23;
T = 300;

doping_na = 1e18;   % p-type (cm^-3)
doping_nd = 1e17;   % n-type (cm^-3)
n_i = 1.5e10;       %(cm^-3)
n_i = n_i*1e6;      % (m^-3)

doping_na = doping_na*1e6;  % (m^-3)
doping_nd = doping_nd*1e6;  % (m^-3)

q = 1.6e-19;

% Defining the region
x = -200:1:200;       % -0.1um to 0.1um
x = x'*nm;              % in m
theta = (x(2) - x(1));

N_A = zeros(size(x));
N_D = zeros(size(x));

N_A(ceil(size(x,1)/2):end) = doping_na;      % p-type region
N_D(1:floor(size(x,1)/2)) = doping_nd;       % n-type region

% Initial Conditions

V = zeros(size(x));
V(ceil(size(x,1)/2):end) = -(k*T/q)*log(doping_na/n_i);
V(1:floor(size(x,1)/2)) = (k*T/q)*log(doping_nd/n_i);
% figure;
% plot(V);
% hold on;
Error = 10; % Arbitrary High value
i = 0;
while Error > 10*eps
    i = i+1;
    d2V_by_dx2=(V(3:end) - 2*V(2:end-1) + V(1:end-2))/theta^2;
    rho = q*(N_D(2:end-1) - N_A(2:end-1) - 2*n_i*sinh(V(2:end-1)/(k*T/q)));
    R = d2V_by_dx2+rho/epsilon;

    Mj = 2/theta^2 + (2*q*n_i/(epsilon*(k*T/q)))*cosh(V(2:end-1,1)/(k*T/q));

    m = size(x,1);
    CM=sparse(1:m-2,1:m-2,Mj,m-2,m-2)...
        +sparse(1:m-2-1,2:m-2,(-1/theta^2)*ones(m-2-1,1),m-2,m-2)+...
        sparse(2:m-2,1:m-2-1,(-1/theta^2)*ones(m-2-1,1),m-2,m-2); 

    DV = CM\R;
    V(2:end-1)=V(2:end-1)+DV;
    Error=norm(DV,2)/sqrt(m);
    % plot(V);
    %plot(rho);
end
figure;
plot(x,V,'LineWidth',1.5);
hold on;
% xlabel('x (m) \rightarrow');
% ylabel('V (volt) \rightarrow');
% title('Potential Profile of the device');
% xlim([x(1),x(end)]);
% line([0 0],[-1 1],'Color','red','LineStyle','--');
% text(-0.5e-7,0,'n-type','HorizontalAlignment','center');
% text(0.5e-7,0,'p-type','HorizontalAlignment','center');

% figure;
% plot(x,-V,'g','LineWidth',2);
% xlabel('x (m) \rightarrow');
% ylabel('Energy (eV) \rightarrow');
% title('Energy Band Diagram');
% xlim([x(1),x(end)]);
% line([0 0],[-1 1],'Color','red','LineStyle','--');
% text(-0.5e-7,0,'n-type','HorizontalAlignment','center');
% text(0.5e-7,0,'p-type','HorizontalAlignment','center');

%% Electron and Hole Densities
n_x = n_i*exp(q*V/(k*T))*1e-6;   % (cm^-3)
h_x = n_i*exp(-q*V/(k*T))*1e-6;   % (cm^-3)
 
% figure;
% plot(x,n_x,'g','LineWidth',1.5);
% hold on;
% plot(x,h_x,'r','LineWidth',1.5);
% xlabel('x (m) \rightarrow');
% ylabel('density   (cm^-^3) \rightarrow');
% title('Electron & Hole Densities');
% xlim([x(1),x(end)]);
% line([0 0],[0 n_x(1)],'Color','blue','LineStyle','--');
% text(-0.5e-7,n_x(1)/2,'n-type','HorizontalAlignment','center');
% text(0.5e-7,n_x(1)/2,'p-type','HorizontalAlignment','center');
% legend('n(x)','h(x)');

% figure;
% semilogy(x,n_x,'g','LineWidth',1.5);
% hold on;
% semilogy(x,h_x,'r','LineWidth',1.5);
% xlabel('x (m) \rightarrow');
% ylabel('ln(density(cm^-^3)) \rightarrow');
% title('Electron & Hole Densities');
% xlim([x(1),x(end)]);
% line([0 0],[0 n_x(1)],'Color','blue','LineStyle','--');
% text(-0.5e-7,exp(log(n_x(1))/2),'n-type','HorizontalAlignment','center');
% text(0.5e-7,exp(log(n_x(1))/2),'p-type','HorizontalAlignment','center');
% legend('n(x)','h(x)');
% hold off;

%% Ideal case calculations
V0 = V(1) - V(end);
xn = sqrt((2*epsilon*V0/q)*doping_na/(doping_nd*(doping_na+doping_nd)));
xp = sqrt((2*epsilon*V0/q)*doping_nd/(doping_na*(doping_na+doping_nd)));
xp_index = knnsearch(x,xp);
xn_index = knnsearch(x,xn);

V_ideal = zeros(size(x));

V_ideal(ceil(size(x,1)/2):xn_index) = (q*doping_nd/epsilon)*(xn-x(ceil(size(x,1)/2):xn_index)/2).*x(ceil(size(x,1)/2):xn_index);
V_ideal(size(x,1)-xp_index:ceil(size(x,1)/2)) = (q*doping_na/(2*epsilon))*(x(size(x,1)-xp_index:ceil(size(x,1)/2))+xp).^2;

V_ideal(xn_index+1:end) = V_ideal(xn_index);
V_ideal(1:ceil(size(x,1)/2)) = V_ideal(1:ceil(size(x,1)/2)) - V_ideal(ceil(size(x,1)/2));
V_ideal = flip(V_ideal,1);

delta = V_ideal(1) - V(1);
V_ideal(:,1) = V_ideal(:,1) - delta;

plot(x,V_ideal,'LineWidth',1.5);
xlabel('x (m) \rightarrow');
ylabel('V (volt) \rightarrow');
title('Potential Profile of the device');
xlim([x(1),x(end)]);
line([0 0],[-1 1],'Color','red','LineStyle','--');
text(-0.5e-7,0,'n-type','HorizontalAlignment','center');
text(0.5e-7,0,'p-type','HorizontalAlignment','center');
legend('Numerical Solution','Depletion Approximation');

close all;clear;
%% NUMERICALLY SOLVING POISSON FOR N-MOSCAP

% Defining constants
epsilon0 = 8.854 * 10^-12;epsilon_si = 11.68;   % S.I. units
epsilon_sio2 = 3.9;
nm = 1e-9;                                      % nano meter -> m
k = 1.38e-23;                                   % S.I. units
q = 1.6e-19;                                    % S.I. units
T = 300;                                        % K

% p-type
doping_na = [1e15,1e16,1e17,1e18]';     %(cm^-3)
doping_na = doping_na*1e6;              %(m^-3)
n_i = 1.5e10;                           %(cm^-3)
n_i = n_i*1e6;                          %(m^-3)

epsilon = epsilon0*epsilon_si;
epsilon_ox=epsilon_sio2*epsilon0;

% Defining the region and parameters inside it
x = 0:1:1000;                      % (nm)
x = x'*nm;                          % (m)
theta = (x(2) - x(1));              % (m)

% Constructing test and result vectors
Vs = [-0.2:0.001:1.5]';
length = size(Vs,1);

figure;
hold on;
for j=1:size(doping_na,1)
    N_A = doping_na(j)*ones(size(x));      % p-type region
    np0 = n_i^2/doping_na(j);
    pp0 = doping_na(j);
    
    V_applied = zeros(length,1);
    charge = zeros(length,1);
    capacitor = zeros(length,1);
    
    % Calculate LFCV for a particular work function
    for i=1:length
        V = zeros(size(x));
        V(1) = Vs(i);   % Defining the Boundary Conditions
        Error = 10;     % Arbitrary High value
        
        % Calculate Potential Profile 
        while Error > 10*eps
            d2V_by_dx2=(V(3:end) - 2*V(2:end-1) + V(1:end-2))/theta^2;
            rho = q*( - N_A(2:end-1) + ( -np0*exp(q*V(2:end-1)/(k*T)) + pp0*exp(-q*V(2:end-1)/(k*T))));
            R = d2V_by_dx2+rho/epsilon;

            Mj = 2/theta^2 + (q/epsilon)*((q/(k*T))*np0*exp(q*V(2:end-1)/(k*T)) + (q/(k*T))*pp0*exp(-q*V(2:end-1)/(k*T)));

            m = size(x,1);
            CM=sparse(1:m-2,1:m-2,Mj,m-2,m-2)...
                +sparse(1:m-2-1,2:m-2,(-1/theta^2)*ones(m-2-1,1),m-2,m-2)+...
                sparse(2:m-2,1:m-2-1,(-1/theta^2)*ones(m-2-1,1),m-2,m-2); 

            DV = CM\R;
            V(2:end-1)=V(2:end-1)+DV;
            Error=norm(DV,2)/sqrt(m);
            %plot(V);
        end

        %rho = q*( - N_A(1:end) + ( -np0*exp(q*V(1:end)/(k*T)) + pp0*exp(-q*V(1:end)/(k*T))));
        E_si = -(V(2)-V(1))/theta;
        E_sio2 = (E_si*epsilon)/epsilon_ox;
        thickness = 10*nm;
        V_applied(i) = Vs(i) + E_sio2*thickness;
        charge(i) = -E_sio2*epsilon_ox;
    end
    
    charge=charge*1e-4;
    
    % Calculating Capacitances
    capacitor(1:end-1) = (charge(2:end)-charge(1:end-1))./(V_applied(2:end) - V_applied(1:end-1));
    capacitor(end) = capacitor(end-1);
    
    plot(V_applied,-capacitor,'LineWidth',1.5);
end

%% Configuring the LFCV Curve
title('LFCV Curve');
xlabel('V_{applied} \rightarrow');
ylabel('C (F/cm^2) \rightarrow');
xlim([-2 9]);
%ylim([0 1.2]);
legend('N_A = 1e15','N_A = 1e16','N_A = 1e17','N_A = 1e18');
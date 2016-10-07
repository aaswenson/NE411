%% 1D Fuel Temperature Analysis

% Alexander Swenson

% NE 411 HW 3

clear; clc; clf;

%% set up model 

void = [0,.1,.14];

for j = 1:length(void)
r_void = void(j);
T_s = 600;                      % clad surface temperature [K]
r_fuel = 0.5;                   % fuel radius [cm]
% r_void = .14;                   % void caused by sintering [cm]
frac_void = (r_fuel^2 - r_void^2)/(r_fuel^2);  % percent of fuel taken by void
gen_multiplier = 1 + frac_void; % increase generation term to account for sintering
g = 145*0.82*gen_multiplier;    % generation in fuel [W/cc]
thick_gap = 0.0089;             % thickness gap [cm]
thick_clad = 0.066;             % thickness clad [cm]
r_gap = r_fuel + thick_gap;
r_clad = r_gap + thick_clad;

kf = .036;                      % fuel conductivity [W/cm-K]
kg = 0.003;                     % gap conductivity [W/cm-K]
kc = .13;                       % clad conductivity [W/cm-K]

fuel_nodes = 100; %nodes(j);          % nodes for each region
gap_nodes =  100; %nodes(j);
clad_nodes = 100;  %nodes(j);
nodes = 100;
total_nodes = fuel_nodes + gap_nodes + clad_nodes;

first_gap_node = fuel_nodes ;
first_clad_node = first_gap_node + gap_nodes;


bound_1 = fuel_nodes;
bound_2 = fuel_nodes + gap_nodes;

dr_fuel = (r_fuel-r_void)/fuel_nodes;
dr_gap = (r_gap - r_fuel)/gap_nodes;
dr_clad = (r_clad - r_gap)/clad_nodes;

r = zeros(total_nodes,1);


% define node points throughout the model

for i = 1:fuel_nodes+1
    r(i) = dr_fuel*(i-1) + r_void;
end

for i = first_gap_node:first_clad_node - 1;
    r(i) = dr_gap*(i-first_gap_node) + r_fuel;
end

for i = first_clad_node : total_nodes 
    r(i) = dr_clad*(i-first_clad_node) + r_gap;
end


%% set up node equations

node_coef = zeros(total_nodes,total_nodes);
b_nodes = zeros(total_nodes,1);

% adiabatic boundary condition at fuel centerline 

node_coef(1,1) = -1 - r(1)/dr_fuel;
node_coef(1,2) = 1 + r(1)/dr_fuel;
b_nodes(1) = -(g*r(1)*dr_fuel)/kf;

for i = 2:bound_1 - 1
    node_coef(i,i) = -1 - (2*r(i))/dr_fuel;
    node_coef(i,i+1)= 1+(r(i)/dr_fuel);
    node_coef(i,i-1) =  r(i)/dr_fuel;
    b_nodes(i) = -(g*r(i)*dr_fuel)/kf;
end 

% couple gap and fuel with two boundary conditions

% BC 1 (use conservation of heat flux to solve for T(n) in fuel


node_coef(bound_1,bound_1-1) = (r(bound_1))/(r(bound_1)-dr_fuel);         % T_n-1 coefficient
node_coef(bound_1,bound_1) = -((r(bound_1)/(r(bound_1)-dr_fuel))+((kg*dr_fuel*(r(bound_1)+dr_gap)*r(bound_1))/(kf*dr_gap*(r(bound_1)-dr_fuel)^2)));
node_coef(bound_1,bound_1+1) = (kg*dr_fuel*(r(bound_1)+dr_gap)*r(bound_1))/(kf*dr_gap*(r(bound_1)-dr_fuel)^2);
b_nodes(bound_1) = 0; 

for i = bound_1 + 1: bound_2 - 1
    
    node_coef(i,i) = -1 - (2*r(i))/dr_gap;
    node_coef(i,i+1)=  1+r(i)/dr_gap;
    node_coef(i,i-1) =  r(i)/dr_gap;
    b_nodes(i)=0;
end


% couple gap and fuel with two boundary conditions
% BC 1 (use conservation of heat flux to solve for T(n) in fuel


node_coef(bound_2,bound_2 -1) = r(bound_2)/(r(bound_2)-dr_gap);       % T_n-1 coefficient
node_coef(bound_2,bound_2) = -((r(bound_2)/(r(bound_2)-dr_gap)) + ((kc*dr_gap*(r(bound_2)+dr_clad)*r(bound_2))/(kg*dr_clad*(r(bound_2)-dr_gap)^2)));
node_coef(bound_2,bound_2+1) = (kc*dr_gap*(r(bound_2)+dr_clad)*r(bound_2))/(kg*dr_clad*(r(bound_2)-dr_gap)^2);

% clad nodes

for i = bound_2 + 1: total_nodes - 1
    node_coef(i,i) = -1 - (2*r(i))/dr_clad;
    node_coef(i,i+1)= 1+r(i)/dr_clad;
    node_coef(i,i-1) = r(i)/dr_clad;
end


% Surface Node
node_coef(total_nodes,total_nodes-1) = r(total_nodes)/dr_clad;
node_coef(total_nodes,total_nodes) = -1 - ((2*r(total_nodes))/dr_clad);
b_nodes(total_nodes) = -T_s*(1+r(total_nodes)/dr_clad);

%% Solve for Temperature Distribution
T = node_coef\b_nodes;
[V,D] = eig(node_coef);

%% Analytical 
T1 =  -.5*(((log(r_gap/r_clad)*g*kg*r_fuel^2)+r_fuel^2*g*log(r_fuel/r_gap)*kc-2*T_s*kc*kg)/(kc*kg));
T2 = -.5*((log(r_gap/r_clad)*r_fuel^2*g-2*kc*T_s)/(kc));
T_a = zeros(total_nodes,1);
for i = 1:total_nodes
    if i <= bound_1 
        T_a(i) = ((r_fuel^2-r(i)^2)*g)/(4*kf) + T1;
    elseif i > bound_1 && i <= bound_2
        T_a(i) = -((T2-T1)*log(r(i)))/(log(r_fuel/r_gap)) + (T2*log(r_fuel)-T1*log(r_gap))/log(r_fuel/r_gap);
    elseif i > bound_2
        T_a(i) = -((T_s-T2)*log(r(i)))/(log(r_gap/r_clad)) + (T_s*log(r_gap)-T2*log(r_clad))/log(r_gap/r_clad);
    end
end
%% Numerical Flux
q = zeros(total_nodes,1);
q(1) = 0;
 for i = 2:bound_1-1
     q(i) = -kf*(T(i+1) - T(i-1))/(2*dr_fuel);
 end
 
 for i = bound_1+1:bound_2-1
     q(i) = -kg*(T(i+1) - T(i-1))/(2*dr_gap);
 end
 
 for i = bound_2+1:total_nodes-1
     q(i) = -kc*(T(i+1) - T(i-1))/(2*dr_clad);
 end
 
 q(bound_1) = (q(bound_1-1)+q(bound_1+1))/2;
 q(bound_2) = (q(bound_2-1)+q(bound_2+1))/2;
 q(total_nodes) = q(total_nodes-1);
 
 %% Analytical Heat Flux
 q_a = zeros(total_nodes,1);
 
 for i = 2:total_nodes
    if i <= bound_1 
        q_a(i) = .5*g*r(i);
    elseif i > bound_1 && i <= bound_2
        q_a(i) = (kg*(T2-T1))/(r(i)*log(r_fuel/r_gap));
    elseif i > bound_2
        q_a(i) = (kc*(T_s-T2))/(r(i)*log(r_gap/r_clad));
    end
 end      

%% Error 
for i= 1:total_nodes
Error_T(i) = (abs(T(i) - T_a(i))/T_a(i))*100;
Error_q(i) = (abs(q(i) - q_a(i))/q_a(i))*100;
end
r_sv{j} = r;
T_sv{j} = T;
q_sv{j} = q;
Err_T_sv{j} = Error_T;
Err_q_sv{j} = Error_q;

end
%% Plotting

colorVec = hsv(length(void));
%Temp

figure(1)
hold on;
for i=1:length(void)
plot(cell2mat(r_sv(i)),cell2mat(T_sv(i)),'color', colorVec(i,:));
legend_info{i} = num2str(void(i));
end
plot(cell2mat(r_sv(length(void))),T_a,'-o')
legend(legend_info,'Analytical Solution')
title('1D Temperature Distribution')
xlabel('radius [cm]')
ylabel('temperature [K]')
hold off;

% Flux

figure(2)
hold on;
for i=1:length(void)
plot(cell2mat(r_sv(i)),cell2mat((q_sv(i))),'color', colorVec(i,:));
legend_info{i} = num2str(void(i));
end
plot(cell2mat(r_sv(length(void))),q_a,'-o')
legend(legend_info,'Analytical Solution')
title('1D Flux Distribution')
xlabel('radius [cm]')
ylabel('Flux [W/cm^2]')
hold off;

% Temperature Error

figure(3)
hold on;
for i=1:length(void)
plot(cell2mat(r_sv(i)),cell2mat((Err_T_sv(i))),'color', colorVec(i,:));
legend_info{i} = num2str(void(i));
end
legend(legend_info)
title('Temperature Profile Error')
xlabel('radius [cm]')
ylabel('% error')
hold off;

% Flux Error 

figure(4)
hold on;
for i=1:length(void)
plot(cell2mat(r_sv(i)),cell2mat((Err_q_sv(i))),'color', colorVec(i,:));
legend_info{i} = num2str(void(i));
end
legend(legend_info)
title('Flux Profile Error')
xlabel('radius [cm]')
ylabel('% error')
hold off;










clear all
clc
tic
%% Notes

% All of the code should be set up ready, just hit run to gain the plots and information. 

% This is the code of Gareth Jones linked with the DE4 Optimisation module. The code is
% split out into separate phases below to align with the work done in the
% report for clarity.

% x(1) = L = Thichness (m)
% x(2) = K = Material thermal conductivity (W/m^2.K) 
% x(3) = A = Surface area (m^2)
% T(1) = T(in) = Temperature of the inside of the radiator (K)
% T(2) = T(out) = Temperature of the outside of the radiator (K)
% H(1) = H(in) = Heat transfer coefficient for water (W/m^2.K)
% H(2) = H(out) = Heat transfer coefficient for air (W'm^2,K)

%Inputs
T1 = 343.15;
T2 = 291.15;
H1 = 500;
H2 = 10;

%% Initial Discovery
A1 = 0.25;

[L, K] = meshgrid(0:0.01:0.1,0:10:500); 

q = -(((T1 - T2)./((1./H1)+(L./K)+(1./H2))).*A1); %Defyning the function

figure;
surf(L,K,q) %Plotting on a 3D Graph

xlabel('Thichkness of Panel (m)');
ylabel('Thermal Conductivity of Material (W/m.K)');
zlabel('Heat Flux (W)');
set(gcf,'color','white')
title('Heat flux, against thickness and thermal conductivity') %Graph styling

%% Interior Point Function

fun = @(x)-(((T1 - T2)/((1\H1) + (x(1)/x(2)) + (1/H2)))*x(3)); %Heat Flux Function
x0 = [0.0015, 120, 1.0]; %Starting position
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.0003, 20, 0.12]; %Lower bounds
ub = [0.0025, 250, 1.6]; %Upper bounds


%SQP
options1 = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunEvals',1000);
[x, fval, exitflag, output, lambda] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@cons,options1)

toc

disp(['Initial Objective: ' num2str(fun(x0))])

disp(table(x(1),x(2),x(3), 'VariableNames',{'Thickness','Material Thermal Conductivity','Surface Area'}))
disp(['Final Objective SQP: ' num2str(fun(x))])

%% Global Search

tic
rng default

opts = optimoptions(@fmincon, 'Algorithm', 'sqp');

problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',@cons,'lb',lb,'ub',ub);
gs = GlobalSearch;
x = run(gs,problem);

disp(table(x(1),x(2),x(3), 'VariableNames',{'Thickness', 'Material Thermal Conductivity', 'Surface Area'}))

disp(['Global Search:' num2str(fun(x))])
toc

%% Non Linear Constraints

function [c, ceq] = cons(x)
T1 = 343.15;
T2 = 291.15;
H1 = 500;
H2 = 10;

ceq = [];
c1 = (((T1 - T2)/((1\H1) + (x(1)/x(2)) + (1/H2)))*x(3));
c = [c1];
end
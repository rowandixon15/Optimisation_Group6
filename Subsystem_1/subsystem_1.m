%% Function finding (based on ASHRAE data)
clear;clc;close all

W = 5;          %Width of wall
H = 10;         %Height of wall
Lw = 0.4;       %Max. thickness of wall
T = 24*60*60;   %Time period of temperature fluctuations (24-hour day)
N = 3;          %Number of layers

data = readtable('Materials Database 2.xlsx');
data.Properties.VariableNames = {'Material' 'L' 'k' 'rho' 'c'};
data.L = data.L./1000; %Convert mm to m

W = 5;          %Width of wall
H = 10;         %Height of wall
Lw = 0.4;       %Max. thickness of wall
T = 24*60*60;   %Time period of temperature fluctuations (24-hour day)

omega = 1/T;    %Angular frequency

a = data.k./(data.rho.*data.c); %Thermal diffusivity
lambda = 2*pi*sqrt(2.*a./omega);%Thermal wavelength

R = data.L./data.k; %Thermal resistance

B = R/(2*pi*(1+1j)*(data.L./lambda))*sinh(2*pi*(1+1j)*(data.L./lambda));
A = cosh(2*pi*(1+1j)*(data.L./lambda));

Z_in = B./(A-1); %Impedance

Cth_in = -1./(omega*imag(Z_in));   %Thermal capacitance

data = addvars(data,R,'After','c');
data = addvars(data,Cth_in,'After','c');

syms X1 X2 X3 X4 X5

%From the data find a single equation for thermal resistance and thermal capacitance

R_eq = X1/X2;  %Polyfit was awful for this, so did it manually

Cth_eqstruct = polyfitn([data.L,data.k,data.rho,data.c],data.Cth_in,2);
Cth_eq = polyn2sym(Cth_eqstruct);

%Lets create another dataset with only the material properties, might come in handy
datamaterial = table(data.Material,data.k,data.rho,data.c,data.R./data.L,data.Cth_in./data.L);
datamaterial.Properties.VariableNames = {'Material' 'k' 'rho' 'c' 'R' 'Cth_in'};

%And remove any repeated materials from it, obviously
index = 2;
while index < length(datamaterial.Material)
    if round(datamaterial.Cth_in(index)) == round(datamaterial.Cth_in(index-1))
        datamaterial(index,:) = [];
    else
        index = index + 1;
    end
end

%% Function checking (based on eyeballing)

R_check = zeros(length(data.L),1);
Cth_check = zeros(length(data.L),1);

%Substitute in variables for each material into calculated equation
for i = 1:length(R_check)
    R_check(i) = subs(R_eq,[X1 X2],[data.L(i) data.k(i)]);
    Cth_check(i) = subs(Cth_eq,[X1 X2 X3 X4],[data.L(i) data.k(i) data.rho(i) data.c(i)]);
end

%Add calculated variables to the data table
data = addvars(data,R_check,'After','R');
data = addvars(data,Cth_check,'After','Cth_in');

%Bar chart for thermal resistance
figure;
bar(categorical(data.Material),[data.R data.R_check]);
legend('Original', 'Calculated');
ylabel('Thermal Resistance  /  m^{2} K W^{-1}');
set(gcf,'color','w');
title('Evaluating Thermal Resistance Function');

%Bar chart for thermal capacitance
figure;
bar(categorical(data.Material),[data.Cth_in data.Cth_check]);
legend('Original', 'Calculated');
ylabel('Thermal Capacitance  /  J m^{-2} K^{-1}');
set(gcf,'color','w');
title('Evaluating Thermal Capacitance Function');

%Scatter chart for both
figure;
scatter(data.R,data.Cth_in,'o');
hold on
scatter(data.R_check,data.Cth_check,'o');
hold off
legend('Original', 'Calculated');
xlabel('Thermal Resistance  /  m^{2} K W^{-1}');
ylabel('Thermal Capacitance  /  J m^{-2} K^{-1}');
% text(data.R,data.Cth_in,"  " + data.Material);
set(gcf,'color','w');
title('Evaluating both Thermal Resistance and Capacitance');

%% Optimise single materials based on normalisation

%Normalise data for each objective
R_thick = data.R./data.L;
Cth_thick = data.Cth_in.*data.L;
R_norm = (R_thick-min(R_thick)) / (max(R_thick)-min(R_thick));
Cth_norm = (Cth_thick-min(Cth_thick)) / (max(Cth_thick)-min(Cth_thick));

%Combined objective function
alpha = 0.5;    %Weighting toward thermal resistance
Normalised = (alpha*R_norm + (1-alpha)*Cth_norm);

%Add normalised variable to the data table
data = addvars(data,Normalised,'After','R_check');

figure;
b = bar(categorical(data.Material),data.Normalised)
ylabel('Normalised Function');
set(gcf,'color','w');
title('Combined and normalised function');
labels1 = string(b(1).YData);
text(labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')

%% Optimise for thermal resistance

%Reformat equation to be optimisable
x = [X1,X2];
funR = matlabFunction(-R_eq,'Vars',{x});

%constraints
x0 = [0 min(data.k)];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [min(data.L) min(data.k)];
ub = [0.4 max(data.k)];
nonlcon = @constraints;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

tic = now;
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(funR, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
% [x,fval,exitflag,output,population,scores] = ga(funR, 2, A, b, Aeq, beq, lb, ub, nonlcon,options);
%Max. L and min. k (obviously)
toc = now;
%% Optimise for thermal capacitance

%Reformat equation to be optimisable
x = [X1,X2,X3,X4];
funCth = matlabFunction(-Cth_eq,'Vars',{x}); 

%constraints
x0 = [0, min(data.k), min(data.rho), min(data.c)];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [min(data.L), min(data.k), min(data.rho), min(data.c)];
ub = [0.4, max(data.k), max(data.rho), max(data.c)];
nonlcon = @constraints;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

tic = now;
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(funCth, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
% [x,fval,exitflag,output,population,scores] = ga(funCth, 4, A, b, Aeq, beq, lb, ub, nonlcon, options);
%Max. L, max. k, max. rho, middly c
toc = now;
%% Combine functions

%Normalise equations for thermal resistance and capacitance
R_eqnorm = (R_eq - min(data.R))/(max(data.R)-min(data.R));
Cth_eqnorm = (Cth_eq - min(data.Cth_in))/(max(data.Cth_in)-min(data.Cth_in));

%Combine functions
x = [X1,X2,X3,X4];
funNorm = matlabFunction(-0.5*R_eqnorm + -0.5*Cth_eqnorm,'Vars',{x});

%constraints
x0 = [0, min(data.k), min(data.rho), min(data.c)];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [min(data.L), min(data.k), min(data.rho), min(data.c)];
ub = [0.4, max(data.k), max(data.rho), max(data.c)];
nonlcon = @constraints;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

tic = now;
% [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(funNorm, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
[x,fval,exitflag,output,population,scores] = ga(funNorm, 4, A, b, Aeq, beq, lb, ub, nonlcon, options);
%Max. L, max. k, max. rho, middly c
toc = now;
%% Optimise both at once with a weighted sum and Pareto frontier

%constraints
x0 = [0, min(data.k), min(data.rho), min(data.c)];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0, min(data.k), min(data.rho), min(data.c)];
ub = [0.4, max(data.k), max(data.rho), max(data.c)];
nonlcon = @constraints;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

n = 10; %number of pareto points
weight = linspace(0,1,n); %weights on each function

%Pareto optimal results
Rvalues = zeros(1,n);
Cthvalues = zeros(1,n);

%Normalise equations for thermal resistance and capacitance
R_eqnorm = (R_eq - min(data.R))/(max(data.R)-min(data.R));
Cth_eqnorm = (Cth_eq - min(data.Cth_in))/(max(data.Cth_in)-min(data.Cth_in));

% R_eqnorm = (R_eq - mean(data.R))/(std(data.R));
% Cth_eqnorm = (Cth_eq - mean(data.Cth_in))/(std(data.Cth_in));

%Optimise combined function at different weightings 
for i = 1:n
    x = [X1,X2,X3,X4];
    %Combine functions at new weightings
    funNorm = matlabFunction(-(weight(i)*R_eqnorm + (1-weight(i))*Cth_eqnorm),'Vars',{x});
    %Optimise combined function
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(funNorm, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    %[x,fval,exitflag,output,lambda,grad,hessian] = ga(funNorm, 4, A, b, Aeq, beq, lb, ub, nonlcon, options);
    %Plug optimal values into original equations and save
    Rvalues(i) = subs(R_eq, [X1 X2], [x(1) x(2)]);
    Cthvalues(i) = subs(Cth_eq,[X1 X2 X3 X4], [x(1) x(2) x(3) x(4)]);
end

figure;
plot(Rvalues, Cthvalues);
hold on
scatter(datamaterial.R*0.4,datamaterial.Cth_in*0.4,'o');
hold off
% text(datamaterial.R*0.4,datamaterial.Cth_in*0.4,"  " + datamaterial.Material);
legend('Pareto set','Real materials at 0.4 m');
xlabel('Thermal Resistance  /  m^{2} K W^{-1}');
ylabel('Thermal Capacitance  /  J m^{-2} K^{-1}');
xlim([0 max(Rvalues)*1.1]);
ylim([0 max(Cthvalues)*1.1]);
set(gcf,'color','w');
title('Theoretical pareto front compared to real materials');

%% Combinations of best thermal capacitance and best thermal insulation

%Best Cth material
[m i] = max(datamaterial.Cth_in);
Best_Cth = datamaterial(i,:);
%Best R material
[m i] = max(datamaterial.R);
Best_R = datamaterial(i,:);

n = 30; %Number of pareto points
weight = linspace(0,1,n);
weight2 = linspace(0,1,n);

Scores2 = zeros(2,n);
for i = 1:n
    Scores2(1,i) = (weight(i)*Best_Cth.R      + (1-weight(i))*Best_R.R)*0.4;      %Thermal resistance
    Scores2(2,i) = (weight(i)*Best_Cth.Cth_in + (1-weight(i))*Best_R.Cth_in)*0.4; %Thermal capacitance   
end

figure;
plot(Rvalues, Cthvalues);
hold on
scatter(datamaterial.R*0.4,datamaterial.Cth_in*0.4,'o');
hold on
scatter(Scores2(1,:),Scores2(2,:),'x');
hold off
% text(datamaterial.R*0.4,datamaterial.Cth_in*0.4,"  " + datamaterial.Material);
legend('Pareto set','Real materials at 0.4 m','Combined Materials');
xlabel('Thermal Resistance  /  m^{2} K W^{-1}');
ylabel('Thermal Capacitance  /  J m^{-2} K^{-1}');
xlim([0 max(Rvalues)*1.1]);
ylim([0 max(Cthvalues)*1.1]);
set(gcf,'color','w');
title('Combination of best thermal resistance and capacitance');

%% Combinations of best thermal capacitance and best thermal insulation and overoptimal

%Best Cth material
[m i] = max(datamaterial.Cth_in);
Best_Cth = datamaterial(i,:);
%Best R material
[m i] = max(datamaterial.R);
Best_R = datamaterial(i,:);
%Best both
Best_Both = datamaterial(5,:); %Wood siding

n = 30; %Number of pareto points
weight = linspace(0,1,n);
weight2 = linspace(0,1,n);

Scores3 = zeros(2,n*n);
for i = 1:n
    for j = 1:n
        Scores3(1,(i-1)*n+j) =(weight(i)*Best_Cth.R      + (1-weight(i))*weight2(j)*Best_R.R       + (1-weight(i))*(1-weight2(j))*Best_Both.R)*0.4;     %Thermal resistance
        Scores3(2,(i-1)*n+j) = (weight(i)*Best_Cth.Cth_in + (1-weight(i))*weight2(j)*Best_R.Cth_in + (1-weight(i))*(1-weight2(j))*Best_Both.Cth_in)*0.4;%Thermal capacitance   
    end
end

figure;
plot(Rvalues, Cthvalues);
hold on
scatter(datamaterial.R*0.4,datamaterial.Cth_in*0.4,'o');
hold on
scatter(Scores3(1,:),Scores3(2,:),'x');
hold off
text(datamaterial.R*0.4,datamaterial.Cth_in*0.4,"  " + datamaterial.Material);
legend('Pareto set','Real materials at 0.4 m','Combined Materials');
xlabel('Thermal Resistance  /  m^{2} K W^{-1}');
ylabel('Thermal Capacitance  /  J m^{-2} K^{-1}');
xlim([0 max(Rvalues)*1.1]);
ylim([0 max(Cthvalues)*1.1]);
set(gcf,'color','w');
title('Combination of materials including the most optimal');

%% nonlcon Function

function [c,ceq]=constraints(x)
    c = [];
    ceq = [];
end

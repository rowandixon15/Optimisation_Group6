clear all
tic
%%% add parameter inputs and obj into a function

% x(1) = Current (Amps)
% x(2) = Length of element (meters)
% x(3) = radius of element (meters)

rb=0.3 % rb= radius of boiler

fun = @(x)(2016000/((((x(1)^2)*x(2)*(0.0000018))/((x(3)^2)*pi))+30000)) %0.0000018 is the resistivity of titanium alloy (selected from CES EduPack)
x0 = [14,1,0.01]; %initial values to minimise from 
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0.001 0.0005]; %lower bounds of varibles (see x1-3 at top)
ub = [45,10,rb];   %upper bounds of varibles (see x1-3 at top)(rb= radius of boiler)

%---------------------N.B.--------------
%Please comment/uncomment out each of the following 3 algorithms to run
%each one.

%---------------------------------- SQP Algorithm ---------------------------------

options1 = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'MaxFunEvals',1000);
[x, fval, exitflag, output, lambda] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@cons,options1)

toc

disp(['Initial Objective: ' num2str(fun(x0))])
disp(table(x(1),x(2),x(3), 'VariableNames',{'current', 'length', 'radius'}))
disp(['Final Objective SQP: ' num2str(fun(x))])

%---------------------------------- Interior Point Algorithm ---------------------------------
% 
% options2 = optimoptions('fmincon','Algorithm', 'interior-point', 'MaxFunEvals',1000);
% 
% [x,fval,ef,output,lambda]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@cons, options2);
% toc
% disp(['Initial Objective: ' num2str(fun(x0))])
% disp(table(x(1),x(2),x(3), 'VariableNames',{'Current', 'Length', 'Radius'}))
% disp(['Final Objective Interior-Point: ' num2str(fun(x))])

%---------------------------------------- Active Set ---------------------------------

% options3 = optimoptions('fmincon','Algorithm', 'active-set', 'MaxFunEvals',1000);
% 
% [x,fval,ef,output,lambda]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@cons, options3);
% toc
% disp(['Initial Objective: ' num2str(fun(x0))])
% disp(table(x(1),x(2),x(3), 'VariableNames',{'Current', 'Length', 'Radius'}))
% disp(['Final Objective Active Set: ' num2str(fun(x))])


%---------------------------------------- Global Search ---------------------------------

rng default

opts = optimoptions(@fmincon,'Algorithm','sqp');

problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',@cons,'lb',lb,'ub',ub);
gs = GlobalSearch;
x = run(gs,problem);

disp(table(x(1),x(2),x(3), 'VariableNames',{'current', 'length', 'radius'}))

disp(['Global Search: ' num2str(fun(x))])


%% Non Linear Constraints

function [c,ceq] = cons(x)
%t = 1:45;
ceq= [];
c1 =  2016000 - ((x(1)^2)*0.0000018*x(2)/(pi*x(3)^2))*45 -(30000*45)  ; %2016000 is the energy needed to heat 12L of water by 40C, 30000 is the power of the gas burner

c = [c1];
end
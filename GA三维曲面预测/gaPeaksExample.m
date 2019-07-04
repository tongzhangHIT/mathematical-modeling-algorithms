%% |peaks| Minimization with Genetic Algorithm
%% Objective Function
% find the minimum of the |peaks| function
clc % ����
clear all; % ɾ��workplace����
close all; % �ص���ʾͼ�δ���

peak

%% Nonlinear Constraint Function
% Subject to a nonlinear constraint defined by a circular region of radius
% three around the origin
type circularConstraint

%% Define Optimization Problem
problem = createOptimProblem('fmincon',...
                             'objective',@(x) peak(x(:,1),x(:,2)), ...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             'nonlcon',@circularConstraint,...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             'x0',[-1 -1],...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                             'lb',[-3 -3],...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             'ub',[3 3],...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             'options',optimset('OutputFcn',...
                                                @peaksPlotIterates))
                             
%% Run the solver |fmincon| from the inital point
% We can see the solution is not the global minimum
%[x,f] = fmincon(problem)                                                    

%% Use Genetic Algorithm to Find the Global Minimum
% Solve the problem using ga.
problem.solver  = 'ga';
problem.fitnessfcn = problem.objective;
problem.nvars = 2;
problem.A=[1,1];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem.b=1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
problem.Aeq=[1,2];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
problem.beq=1.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
problem.options = gaoptimset('PopInitRange',[-3;3],...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             'OutputFcn',@peaksPlotIterates,...
                             'Display','iter',[1,2])

[x,f] = ga(problem)
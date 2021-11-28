function [ Jirr, Jflo ] = evaluate_objective(x, M, V)
%
% function f = evaluate_objective(x, M, V)
%
% Function to evaluate the objective functions for the given input vector x.%
% x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables. 
%
% This functions is basically written by the user who defines his/her own
% objective function. Make sure that the M and V matches your initial user
% input.
%

x = x(1:V) ;
x = x(:)   ;

% --------------------------------------
% insert here your function:
% global variable to pass inside extra inputs
global opt_inputs ;
n = opt_inputs.n ;
h_init = opt_inputs.h_init ;
param = opt_inputs.param;
h_flo = opt_inputs.h_flo; 
w = param.reg.w ;
% 1) policy param
param.reg.h1 = x(1) ;
param.reg.h2 = x(2) ;
param.reg.h3 = x(3) ;
param.reg.m1 = x(4) ;


% 2) run simulation
[s_reg, h_reg, r_reg] = simulate_reg_lake( n, h_init, param ) ;
% 3) compute 2 objs
h_reg = h_reg(2:end);
r_reg = r_reg(2:end);
% I1: daily average squared deficit
dr = max( w-r_reg, 0 ) ;
Jirr = mean( dr.^2 ) ;
% IF1 = mean flooded indication in the city of Tudela
Ny = length(h_reg)/365;
Jflo = sum( h_reg>h_flo )/Ny

end




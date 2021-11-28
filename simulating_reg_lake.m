function [s, h, r] = simulating_reg_lake( n, h_init, param ) 
delta = 60*60*24;
H = length(n)-1 ; % 10 years horizon
% initialization of vectors
s = nan( size(n) );
h = nan( size(n) );
r = nan( size(n) );
% initial condition t=1
h(1) = h_init; 
s(1) = h_init * param.nat.S ;
for t=1:H
   % 1)release
   r(t+1) = regulating_release( param , h(t) );
   %r(t+1) = param.nat.beta*( h(t) - param.nat.h0 ).^param.nat.alpha;  % this replaces the previous line in case of simulation of the natural lake
   % 2) mass-balance
   s(t+1) = s(t) + ( n(t+1) - r(t+1) )*delta ;
   % 3) s->h
   h(t+1) = s(t+1)/param.nat.S ;
end

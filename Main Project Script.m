
%% Part 2
%Numerical comparison of the performance attained by few
%project alternatives on different stakeholders’ interests



%% PART 2.1: performance of Alternative-0 (i.e. no dam)

load -ascii Tudela.csv
q = Tudela(:,5); % average daily stream flow[m3/s]
h_init = 0.015 ; % initial condition
u = Tudela(:,4); %average daily precipitation[mm/d]
t = Tudela(:,6); %average daily temperature[c]
w = 65; %  assume 65 m3/s as water demand for irrigation
h_flo = 0.04; % assume 0.04 m as flood threshold

% natural level-discharge relationship:
param.nat.S = 25729870000 ; % [m2]
param.nat.alpha = 20000 ; %[m2/s]
param.nat.h0 = 0 ;   %[m]


%reference of natural lake (since Lake Tudela  was a natural system)
n = q;
n = [ nan; n ] ; % for time convention
[s_nat, h_nat, r_nat] = simulating_nat_lake( n, h_init, param ) ;

figure; 
plot(h_nat); ylabel('level (m)'); legend('natural')
figure;
plot(r_nat); ylabel('release (m3/2)'); legend('natural')

h_nat = h_nat(2:end);
r_nat = r_nat(2:end);

% I1: daily average squared deficit
dn = max( w-r_nat, 0 ) ;
I1_nat = mean( dn.^2 ) ;    %  Since our level value is between 0 and 1, 
                            %  the area indicator becomes negative value,
                            % therefore we dont use area indicator
                            % {S _flo( idx ) = 0.081*h_reg(idx).^3
                            % -0.483*h_reg(idx).^2 +1.506*h_reg(idx)-1.578}
                            
                        
                             

% IF1 = mean flooded surface in the city of Locarno
Ny = length(h_nat)/365;
I2_nat = sum( h_nat>h_flo )/Ny;


%% PART 2.2: Active capacity of the dam

% monthly mean flow
T = 27;  %number of years
qMonth = dailyToMonthly(q, T); % m3/s

% target release = downstream water demand w
w = 100 ; % m3/s

% Sequent Peak Analysis
deltaT = 3600*24*[31 28 31 30 31 30 31 31 30 31 30 31]';
Q = qMonth(:).*repmat(deltaT,27,1) ; % m3/month
W = w*ones(size(Q)).*repmat(deltaT,27,1) ; % m3/month
K = zeros(size(Q));
K(1) = 0;

for t = 1:length(Q)
    K(t+1) = K(t) + W(t) - Q(t) ;
    if K(t+1) < 0
        K(t+1) = 0 ;
    end        
end

figure; plot( K );xlabel('time[month]'); ylabel('deficit[m3/month]');
Kopt = max(K);
% once the capacity is defined, it's necessary to set dam height/surface 
% In the case of Lake Tudela, the surface is 25729.870000 km2
% so we can make the following check:
h_max = Kopt/25729870000;
% this dam height value is in the range of different level(comparison to the Lab case_studies) because of it's big area,  
% suggesting that the storage capacity of Lake Tudela, although designed 
% by nature and not artificially, is suitable for providing water supply 
% to the downstream users






%% PART 2.3.a: simulation of water reservoir under regulation

% regulated level-discharge relationship:
param.reg.w = 65 ;   %[m3/s]
param.reg.h_min = 0 ; %[m]
param.reg.h_max = h_max ; %h_max is the height of active dam [m]
param.reg.h1 = 0.01 ; %[m]
param.reg.h2 = 0.02 ; %[m]
param.reg.m1 = 3000 ;  %[m2/s]
param.reg.m2 = 35000 ; %[m2/s]
param.reg.h3 = 0.005 ; %[m]


% test regulated level-discharge relationship
h_test = [ 0:0.002:0.05 ] ;
r_test = regulating_release( param , h_test );
figure; plot(h_test,r_test, 'o');
xlabel('h_test[m]'); ylabel('r_test[m3/s]');

% simulation of lake dynamics
n = q;
n = [ nan; n ] ; % for time convention
h_init = 0.015 ; % initial condition


[s, h, r] = simulating_lake( n, h_init, param ) ;

figure; plot( h ); xlabel('time[d]'); ylabel('water level[m]');

figure; plot(r) ; xlabel('time[d]'); ylabel('release[m3/s]');

% impacts
h1 = h(2:end);
r1 = r(2:end);
w = param.reg.w ;
h_flo = 0.04 ;

% daily average squared deficit
def = max( w-r1, 0 );
jirr = mean( def.^2 );

% daily average flooded area
Ny = length(h1)/365;
jflo = sum( h>h_flo )/Ny;



%% PART 2.3.b: Performance of few management alternatives defined by different parameterizations of the Standard Operating Policy
 
% evaluating different alternatives in comparison with no-dam alternatives(natural lake)

global opt_inputs;
opt_inputs.n = n ;
opt_inputs.h_init = h_init ;
opt_inputs.param = param ;
opt_inputs.h_flo = h_flo ;

[ Jirr, Jflo ] = evaluating_objective([0.01 0.02 0.005 3000 35000], 2, 5) ;     % existing policy indicator
[ Jirr1, Jflo1 ] = evaluating_objective([0.01 0.033 0 3000 35000], 2, 5)  ;      % expected improving water supply indicator
[ Jirr2, Jflo2 ] = evaluating_objective([0.01 0.011 0 5000 50000], 2, 5)  ; % expected improving flood indicator

figure; plot( I1_nat,I2_nat , 'bo') ; 
hold on; plot(  Jirr, Jflo  , 'yo') ; 
hold on; plot(  Jirr1, Jflo1 , 'ro') ;
hold on; plot(  Jirr2, Jflo2 , 'mo') ;
xlabel('irrigation'); ylabel('flood');
legend('nat','reg', 'reg1', 'reg2');


% the alternative-0(no-dam) has the performance of [ I1_nat(irrigation) =
% 144.21 and I2_nat(flo) = 5.25 ]. We manually use different parameter
% value by changing our decision variable(h1 h2 h3 m1 m2) based on
% different stake_holder interest2(flood, irrigation). The existing policy
% with decision variable of ([h1=0.01,h2=0.02, h3=0.005,
% m1=3000,m2=35000]) has the indicator of [Jirr = 2.37 and Jflo = 10.44]
% ,showing that we have high improvement in irrigation and a little worsening
% of performance of flood. The second policy
% with decision variable of ([h1=0.01,h2 =0.033, h3=0,
% m1=3000,m2=35000]) has the indicator of [Jirr1 = 0 and Jflo1 = 13.296]
% ,showing that we have a high improvement in irrigation(which is Jirr1 = 0; this is very interesting!!!!!????) and a considerable worsening
% of performance of flood(Jflo1 = 13.296) in comparison to
% alternative_0(Jflo = 5.25). The third policy
% with decision variable of ([h1=0.01,h2 =0.011, h3=0.005,
% m1=5000,m2=50000]) has the indicator of [Jirr2 =56.03 and Jflo2 = 5.96]
% ,showing that we have very little worsening performance in flood (which
% is Jflo2 = 5.96)in comparison with alternative_0(J_flo = 5.25) and a high
% improvement of performance in irrigation(Jirr2 =56.03) in comparison with alternative_0(J_irr = 144.21).










clc ; clear all; close all;
% SIR Model parameter
xmin=2.25e-8;
xmax=5.4e-8;
beta=xmin+rand(1,1)*(xmax-xmin); % Randomized, rate of infection  R= (N*beta)/gamma
gamma =0.6; % Rate of recovery 
delta = 1/60; % Rate of immunity loss   1/day

N = 40*10^6; % Total population N = S + I + R , e.g.-->(Malaysia estimated population as of 2020, 4*10^7 --adding also non-citizens such as working immigrants)
I0 = 10; % Initial number of infected
T = 365; % Period of 365 days
dt = 1/4; % Time interval of 6 hours (1/4 of a day)
fprintf('R0 value is %.2f',N*beta/gamma)

tt = 0:dt:T-dt;
% Calculate the model
[S,I,R] = sir_model(beta,gamma,delta,N,I0,T,dt);

% Plot SIR Model
figure
plot(tt,S,'b',tt,I,'r',tt,R,'g','LineWidth',3); grid on;
set(gca,'FontSize', 20)
title("Case of Covid-19 in Malaysia",'FontSize',24)
xlabel('Days','FontSize', 20);ylabel('Number of individuals','FontSize', 20);
legend('Susceptible','Infected','Removed','FontSize', 20);
 
function [S,I,R] = sir_model(beta,gamma,delta,N,I0,T,dt)
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    for tt = 1:(T/dt)-1
        % Equations of the model
        dS = (-beta*I(tt)*S(tt) + delta*R(tt)) * dt;
        dI = (beta*I(tt)*S(tt) - gamma*I(tt)) * dt;
        dR = (gamma*I(tt) - delta*R(tt)) * dt;
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
    end
end

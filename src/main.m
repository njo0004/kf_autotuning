clear all; clc; close all;

addpath(genpath("helper_functions"));

%{

    I will be implementing the paper:

    Weak in the NEES?: Auto-Tuning Kalman Filter Using Bayesian
    Optimization

%}

A = [0 1;0 0];
B = [0;1];
H = [1 0];
Gamma = [0;1];

dt = 0.1;
dt_sim = 0.001;

V = 1;
W = 0.1;

signal_amp = 2;
signal_freq = 0.75;

%% Truth Simulation 

T = 0:dt_sim:20; % <- time steps of simulation
N = 20;  % <- number of simulations
meas_idx = 1;

for i = 1:N

    X(:,1,i) = [0;0] + 0.1*randn(2,1); % <- generate IC for each simulation

    for j = 1:length(T)

        Xdot = A*X(:,j,i) + B*signal_amp*sin(j*dt*signal_freq) + Gamma*V*randn;
        X(:,j+1,i) = X(:,j,i) + Xdot*dt_sim;

        if(~mod(dt,dt_sim*j))

            Y(i,meas_idx) = H*X(:,j,i) + (W/dt)*randn(1);
            
            meas_idx = meas_idx + 1;

        end

    end

    meas_idx = 1;

    figure(1)
    plot(T,X(1,1:length(T),i))
    hold on
    xlabel('Time (s)')
    ylabel('Position (m)')
    title('1D Cart Postion')

    figure(2)
    plot(T,X(2,1:length(T),i))
    hold on
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    title('1D Cart Velocity')

    figure(3)
    plot(Y(i,:))
    hold on
    xlabel("Measurement IDX")
    ylabel("1D Position Measurement (m)")
    title("Position Measurement")

end
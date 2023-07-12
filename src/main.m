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
dt_sim = 0.1;

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

        u(i,j) = signal_amp*sin(j*dt*signal_freq);

        Xdot = A*X(:,j,i) + B*u(i,j) + Gamma*V*randn;
        X(:,j+1,i) = X(:,j,i) + Xdot*dt_sim;

        if(~mod(dt_sim*j,dt))
            
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
    grid on

    figure(2)
    plot(T,X(2,1:length(T),i))
    hold on
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    title('1D Cart Velocity')
    grid on

    figure(3)
    plot(Y(i,:))
    hold on
    xlabel("Measurement IDX")
    ylabel("1D Position Measurement (m)")
    title("Position Measurement")
    grid on

end


%% Running Kalman Filters

F = [1 dt;0 1];
G = [0.5*dt^2;dt];
Q = VanLoanDiscretization(V,B,A,dt);
R = W/dt;

X_hat(:,1,:) = zeros(2,1,i);

for i = 1:N
    
    P(:,:,1,i) = 5*eye(2); % OOF

    for j = 1:20/dt

        X_hat(:,j+1,i) = F*X_hat(:,j,i) + G*signal_amp*sin(j*dt*signal_freq);
        
        delX(:,j,i) = X(:,j,i) - X_hat(:,j,i);

        P(:,:,j+1,i) = F*P(:,:,j,i)*F' + Q;
        
        NEES(i,j) = delX(:,j,i)'*inv(P(:,:,j,i))*delX(:,j,i);

        S(:,:,j,i) = (H*P(:,:,j,i)*H' + R);

        K = P(:,:,j+1,i)*H'/(S(:,:,j,i));

        delY(i,j) = Y(i,j) - H*X_hat(:,j+1,i);

        NIS(i,j) = delY(i,j)'/(S(:,:,j,i))*delY(i,j);

        X_hat(:,j+1,i) = X_hat(:,j+1,i) + K*delY(i,j);

        P(:,:,j+1,i) = (eye(2) - K*H)*P(:,:,j+1,i);

    end
    
    figure(3)
    plot(NEES(i,:))
    hold on
    title("NEES")

    figure(4)
    plot(NIS(i,:))
    hold on
    title("NIS")

end

NEES_ave = (1/N).*sum(NEES,1);
NIS_ave  = (1/N).*sum(NIS,1);

figure(3)
plot(NEES_ave)

figure(4)
plot(NIS_ave)

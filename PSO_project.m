clc;
clear all;
close all;
% Lower bound for acceleration and steering rate 

lb = [-0.75 -pi/4];

%Upper bound for acceleration and steering rate 

ub = [0.75 pi/4];
Np = 50;        %Population size
MaxIter = 500;  %Number of iterations
w = 0.9;        %Inertia weight
c1 = 2;         %Acceleration coefficient
c2 = 2;         %Acceleration coefficient
var = 2;        %Control input variables
N = 50;         %Number of time steps
dt = 0.1;       %Time

%%
start_cord = [10.7, 1.95];
dest_cord = [1, -1];

d_ref = ((start_cord(1) - dest_cord(1))^2 + ((start_cord(2) - dest_cord(2))^2))^0.5

x_init = [10.7 1.95 0 -1 -pi/8]; %Initial state1
% x_init = [10.7 1.95 pi/15 -1 -pi/8]; %Initial state2
% x_init = [10.7 1.95 0 -1 -pi/8]; %Initial state3
% x_init = [10.7 1.95 -pi/9 -1 -pi/8]; %Initial state4
L = 2.5;        %length between the front axle and the rear axle

%Initialize a random population within bounds
for i = 1:Np
    X(:,:,i) = repmat(lb, N, 1) + repmat((ub - lb), N, 1).*rand(N, var);
end

%Initialize a random velocity within bounds
for i = 1:Np
    V(:,:,i) = repmat(lb, N, 1) + repmat((ub - lb), N, 1).*rand(N, var);
end

for i = 1:Np
    [~,J(i)] = obj_func_PSO(N, dt, X(:,:,i), L, x_init, d_ref);
end

%Initialize the personal best position
X_pbest = X;
%Initialize the personal best function value
J_pbest = J;
%Determine the global best function value
[J_gbest,index]  = min(J_pbest);
%Determine the best global position
X_gbest = X(:,:,index);

for i = 1:MaxIter
    for j = 1:Np
        %Update the velocity of each particle
        V(:,:,j) = (w*V(:,:,j)) + (c1*rand(N,var).*(X_pbest(:,:,j) - X(:,:,j))) + (c2*rand(N,var).*(X_gbest - X(:,:,j)));
        %Update the position of each particle
        X(:,:,j) = X(:,:,j) + V(:,:,j);
        %Bound the violating variables to their lower bound
        X(:,:,j) = max(X(:,:,j), lb);
        %Bound the violating variables to their upper bound
        X(:,:,j) = min(X(:,:,j), ub);
        %Calculate the function value for each particle
        [~,J(j)] = obj_func_PSO(N, dt, X(:,:,j), L, x_init, d_ref);
        %Check if the current function value if less than the personal best function value
        if (J(j) < J_pbest(j))
            %Update the personal best function value and position
            J_pbest(j) = J(j);
            X_pbest(:,:,j) = X(:,:,j);
            %Check if the personal best function value is less than the global best function value
            if (J_pbest(j) < J_gbest)
                %Update the global best function value and position
                J_gbest = J_pbest(j);
                X_gbest = X_pbest(:,:,j);
            end
        end 
    end
end

J_gbest;
X_gbest;

%% Plotting the near optimal solution from PSO

[xf,~] = obj_func_PSO(N, dt, X_gbest, L, x_init, d_ref);
plot(xf(:,1), xf(:,2), '--b')
xlim([-2 15])
ylim([-5 4])

%% Evaluating the Objective function (Includes penalty functions to tackle constraints)
function [x,f] = obj_func_PSO(N, dt, Xp, L, x_init, d_ref)

    k1 = 1;
    k2 = 2;
    k3 = 2;
    k4 = 2;
    k5 = 2;
    
    for i = 1:N
        if i==1
            x(i,1) = x_init(1) + x_init(4)*cos(x_init(3))*dt;
            x(i,2) = x_init(2) + x_init(4)*sin(x_init(3))*dt;
            x(i,3) = x_init(3) + ((x_init(4)*tan(x_init(5)))/L)*dt;
            x(i,4) = x_init(4) + Xp(i,1)*dt;
            x(i,5) = x_init(5) + Xp(i,2)*dt;
        else
            x(i,1) = x(i-1,1) + x(i-1,4)*cos(x(i-1,3))*dt;
            x(i,2) = x(i-1,2) + x(i-1,4)*sin(x(i-1,3))*dt;
            x(i,3) = x(i-1,3) + ((x(i-1,4)*tan(x(i-1,5)))/L)*dt;
            x(i,4) = x(i-1,4) + Xp(i,1)*dt;
            x(i,5) = x(i-1,5) + Xp(i,2)*dt; 
        end
    end
    
    for z = 1:N-1
        d(z) = ((x(z,1) - x(z+1,1))^2 + (x(z,2) - x(z+1,2))^2)^0.5;
    end
    
    f = k1*abs((d_ref-sum(d))) + k2*(x(N,1)-1)^2 + k3*(x(N,2)+1)^2 + k4*(Xp(N,2))^2 + k5*(Xp(N,1))^2;

end
%% ECE_592 Project  Multi Stage Optimization Framework

%   1. First stage  - Particle Swarm Optimization 
%   2. Second stage - Sequential Quadratic Programming (SQP)

clc;
clear;
close all;
opengl('save', 'software')

%% Calling the Particle Swarm Optimization function to get a near optimal solution

PSO_project  %% Comment it for SQP alone

% Using the near optimal solution from the PSO as an initial guess for SQP for quick convergence 

u_init = X_gbest;

%% Parking lot parameters
% All values are in meters (m)
 
Cl = 4;                    %%% width of the road
Lsw = 2;                   %%% width of the parking space
Lsl = 7.3;                 %%% Length of the parking space

L = 2.5;                   %%% Wheel base
m = 0.7;                   %%% Rear over hang length
n =0.8;                    %%% Front over hang length
b = 0.75;                  %%% Half of trackwidt h  T = 2b =1.7 m


%% Problem setting

%%% States - X = [x,y,theta, v, phi]   

% x         - x coordinate of the reference point
% y         - y coordinate of the reference point
% theta     - orientation of the vehicle's longitudinal axis with the
             % Global reference frame's X axis
% v         - Longitudinal velocity of the vehicle 
% phi       - Steering angle


N = 50;                                 % Time Steps for discretized vehicle model
dt_init = [0.05 0.02];                  % Initial guess for the time step

a_init = -0.75*ones(N,1);               %%% initial acceleration input
omega_init = pi/8 *ones(N,1);           %%% initial steering rate input

x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state1
% x_init = [10.7 1.95 pi/15 -1 -pi/8];  %%% Initial state2
% x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state3
% x_init = [10.7 1.95 -pi/9 -1 -pi/8];  %%% Initial state4


%% Uncomment the following line to check the solution of a Stand alone SQP optimization 

u_init = [a_init omega_init]; 

%% Visualizing the initial position of the car

Ax_0 = x_init(1) + (L+n)*cos(x_init(3))-b*sin(x_init(3));
Ay_0 = x_init(2) + (L+n)*sin(x_init(3))+b*cos(x_init(3));

Bx_0 = x_init(1) + (L+n)*cos(x_init(3))+b*sin(x_init(3));
By_0 = x_init(2) + (L+n)*sin(x_init(3))-b*cos(x_init(3));

Cx_0 = x_init(1) - m*cos(x_init(3))+b*sin(x_init(3));
Cy_0 = x_init(2) - m*sin(x_init(3))-b*cos(x_init(3));

Dx_0 = x_init(1) - m*cos(x_init(3))-b*sin(x_init(3));
Dy_0 = x_init(2) - m*sin(x_init(3))+b*cos(x_init(3)); 


%%% Marking the reference point on the vehicle

plot(x_init(1),x_init(2),'*')

%%%  Limiting the x axis and y axis to better represent the parking lot
%%%  under consideration

xaxis = [-5 15];
yaxis = [-4 4.5];

xlim(xaxis)  
ylim(yaxis)
hold on

%%% Plotting the corner points of the car to locate it 

xCar =[Ax_0 Bx_0 Cx_0 Dx_0 Ax_0];
yCar = [Ay_0 By_0 Cy_0 Dy_0 Ay_0];
plot(xCar,yCar,'--');

%%%                  For plotting the parking space 

%                       (Parking space 1)
% 
plot(0:Lsl,zeros(1,length(0:Lsl)),'--r','LineWidth',2)   
plot(0:Lsl,-Lsw*ones(1,length(0:Lsl)),'--r','LineWidth',2)

%%% for plotting the vertical lines denoting the ends of the parking space| |

plot(zeros(1,length(-Lsw:0)),-Lsw:0,'--r','LineWidth',2)
plot(Lsl*ones(1,length(-Lsw:0)),-Lsw:0,'--r','LineWidth',2)



%%%                     (Parking space 2)
% 
% plot(0:Lsl,zeros(1,length(0:Lsl)),'--r','LineWidth',2)   
% plot(0:Lsl,Lsw*ones(1,length(0:Lsl)),'--r','LineWidth',2)
% 
% %%%for plotting the vertical lines denoting the ends of the parking space| |
% plot(zeros(1,length(0:Lsw)),0:Lsw,'--r','LineWidth',2)
% plot(Lsl*ones(1,length(0:Lsw)),0:Lsw,'--r','LineWidth',2)

% For plotting the horizontal lines --- To denote road

plot(xaxis(1):xaxis(2),Cl*ones(1,length(xaxis(1):xaxis(2))),'--r','LineWidth',1)
plot(xaxis(1):0,0*ones(1,length(xaxis(1):0)),'--r','LineWidth',1)
plot(Lsl:xaxis(2),0*ones(1,length(Lsl:xaxis(2))),'--r','LineWidth',1)

xlabel('X (m)')
ylabel('Y (m)')
%% Calling FMINCON (SGP solver)

% options=optimoptions('fmincon','Display','iter','Algorithm','sqp','maxiterations',1000,'MaxFunctionEvaluations',100000,'OptimalityTolerance',0.0000001);
options=optimoptions('fmincon','Algorithm','sqp','maxiterations',1000,'MaxFunctionEvaluations',100000,'OptimalityTolerance',0.0000001);
[u_opt, J_opt] = fmincon(@obj_func,[u_init;dt_init],[],[],[],[],[],[],@const,options);



%% Plotting results

a_opt = u_opt(1:N,1);
w_opt = u_opt(1:N,2);
dt = u_opt(end,1);
x = zeros(N,5);
  
% Updating the states at each time step using the update equations derived
% by transcribing the vehicle dynamics equations in Continuous time to discrete time
% using Forward Euler approximation

for i = 1:N
    if i==1
        x(i,1) = x_init(1) + x_init(4)*cos(x_init(3))*dt;
        x(i,2) = x_init(2) + x_init(4)*sin(x_init(3))*dt;
        x(i,3) = x_init(3) + ((x_init(4)*tan(x_init(5)))/L)*dt;
        x(i,4) = x_init(4) + a_opt(i)*dt;
        x(i,5) = x_init(5) + w_opt(i)*dt;
    else
        x(i,1) = x(i-1,1) + x(i-1,4)*cos(x(i-1,3))*dt;
        x(i,2) = x(i-1,2) + x(i-1,4)*sin(x(i-1,3))*dt;
        x(i,3) = x(i-1,3) + ((x(i-1,4)*tan(x(i-1,5)))/L)*dt;
        x(i,4) = x(i-1,4) + a_opt(i)*dt;
        x(i,5) = x(i-1,5) + w_opt(i)*dt; 
    end
end
% Pre allocating the values of the corners points of the vehicle at each
% time step

Ax = zeros(N,1);
Bx = zeros(N,1);
Cx = zeros(N,1);
Dx = zeros(N,1);
Ay = zeros(N,1);
By = zeros(N,1);
Cy = zeros(N,1);
Dy = zeros(N,1);
  
 %Finding the corner points of the vehicle from the states using the
 %following Kinematics equations
 
for i =1:N
    Ax(i) = x(i,1) + (L+n)*cos(x(i,3))-b*sin(x(i,3));
    Ay(i) = x(i,2) + (L+n)*sin(x(i,3))+b*cos(x(i,3));

    Bx(i) = x(i,1) + (L+n)*cos(x(i,3))+b*sin(x(i,3));
    By(i) = x(i,2) + (L+n)*sin(x(i,3))-b*cos(x(i,3));

    Cx(i) = x(i,1) - m*cos(x(i,3))+b*sin(x(i,3));
    Cy(i) = x(i,2) - m*sin(x(i,3))-b*cos(x(i,3));

    Dx(i) = x(i,1) - m*cos(x(i,3))-b*sin(x(i,3));
    Dy(i) = x(i,2) - m*sin(x(i,3))+b*cos(x(i,3)); 
end
    
%%% Marking the reference point on the vehicle

figure(2)
plot(x_init(1),x_init(2),'*k')
hold on
plot(x(N,1),x(N,2),'*k')

%%%  Limiting the x axis and y axis to better represent the parking lot
%%%  under consideration

xaxis = [-5 15];  
yaxis = [-4 4.5];

xlim(xaxis)
ylim(yaxis)

hold on

%%% Plotting the corner points of the car in its initial state
xCar =[Ax_0 Bx_0 Cx_0 Dx_0 Ax_0];
yCar = [Ay_0 By_0 Cy_0 Dy_0 Ay_0];
plot(xCar,yCar,'--k');

%%% Plotting the corner points of the car to locate it 

xCar_N =[Ax(N) Bx(N) Cx(N) Dx(N) Ax(N)];
yCar_N = [Ay(N) By(N) Cy(N) Dy(N) Ay(N)];
plot(xCar_N,yCar_N,'--');
hold on


xCar_total =zeros(N/2,5);
yCar_total =zeros(N/2,5);

%%% Plotting the vehicle's profile at each time step 

for i =1: N/2
    xCar_total(i,:) =  [Ax(2*i) Bx(2*i) Cx(2*i) Dx(2*i) Ax(2*i)];
    yCar_total(i,:)=  [Ay(2*i) By(2*i) Cy(2*i) Dy(2*i) Ay(2*i)];
    plot(xCar_total(i,:),yCar_total(i,:),'--b');
    hold on;
end

%%%                  For plotting the parking space 

%                       (Parking space 1)

plot(0:Lsl,zeros(1,length(0:Lsl)),'--r','LineWidth',2)   
plot(0:Lsl,-Lsw*ones(1,length(0:Lsl)),'--r','LineWidth',2)

% % for plotting the vertical lines denoting the ends of the parking space| |

plot(zeros(1,length(-Lsw:0)),-Lsw:0,'--r','LineWidth',2)
plot(Lsl*ones(1,length(-Lsw:0)),-Lsw:0,'--r','LineWidth',2)

x = [x_init;x];

%%%                     (Parking space 2)
% % 
% plot(0:Lsl,zeros(1,length(0:Lsl)),'--r','LineWidth',2)   
% plot(0:Lsl,Lsw*ones(1,length(0:Lsl)),'--r','LineWidth',2)
% 
% % for plotting the vertical lines denoting the ends of the parking space| |
% plot(zeros(1,length(0:Lsw)),0:Lsw,'--r','LineWidth',2)
% plot(Lsl*ones(1,length(0:Lsw)),0:Lsw,'--r','LineWidth',2)

% For plotting the horizontal lines --- To denote road

plot(xaxis(1):xaxis(2),Cl*ones(1,length(xaxis(1):xaxis(2))),'--r','LineWidth',1)
plot(xaxis(1):0,0*ones(1,length(xaxis(1):0)),'--r','LineWidth',1)
plot(Lsl:xaxis(2),0*ones(1,length(Lsl:xaxis(2))),'--r','LineWidth',1)

%%% Plotting the progression of the reference point
plot(x(:,1),x(:,2),'--k')


xlabel('X (m)')
ylabel('Y (m)')

%%% Plotting the states other than x and y
figure

time = linspace(0,N*dt,N+1);
subplot(3,2,1)
plot(time,rad2deg(x(:,3)),'b')
xlim([0,time(end)])
xlabel('Time (s)')
ylabel('Orientation(deg)')

subplot(3,2,2)
plot(time,rad2deg(x(:,5)),'b')
xlim([0,time(end)])
xlabel('Time (s)')
ylabel('Steering angle(deg)')

subplot(3,2,3)
plot(time,x(:,4),'b')
xlim([0,time(end)])
xlabel('Time (s)')
ylabel('Velocity(m/s)')

time = linspace(0,N*dt,N);
subplot(3,2,4)
plot(time,u_opt(1:N,1),'b')
xlim([0,time(end)])
xlabel('Time (s)')
ylabel('Acceleration(m/sec^2)')

subplot(3,2,5)
plot(time,rad2deg(u_opt(1:N,2)),'b')
xlim([0,time(end)])
xlabel('Time (s)')
ylabel('Steering rate(deg/sec)')



%% Objective function -- Here it is the Total time for parking maneuver

function J = obj_func(dec_vars)
    N = 50;
    dt = dec_vars(end,1);
    a_opt = dec_vars(1:N,1);
    w_opt = dec_vars(1:N,2);
    L = 2.5;
    x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state1
%     x_init = [10.7 1.95 pi/15 -1 -pi/8];  %%% Initial state2
%     x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state3
%     x_init = [10.7 1.95 -pi/9 -1 -pi/8];  %%% Initial state4
    x = zeros(N,5);
    
% Updating the states at each time step using the update equations derived
% by transcribing the vehicle dynamics equations in Continuous time to discrete time
% using Forward Euler approximation

    for i = 1:N
        if i==1
            x(i,1) = x_init(1) + x_init(4)*cos(x_init(3))*dt;
            x(i,2) = x_init(2) + x_init(4)*sin(x_init(3))*dt;
            x(i,3) = x_init(3) + ((x_init(4)*tan(x_init(5)))/L)*dt;
            x(i,4) = x_init(4) + a_opt(i)*dt;
            x(i,5) = x_init(5) + w_opt(i)*dt;
        else
            x(i,1) = x(i-1,1) + x(i-1,4)*cos(x(i-1,3))*dt;
            x(i,2) = x(i-1,2) + x(i-1,4)*sin(x(i-1,3))*dt;
            x(i,3) = x(i-1,3) + ((x(i-1,4)*tan(x(i-1,5)))/L)*dt;
            x(i,4) = x(i-1,4) + a_opt(i)*dt;
            x(i,5) = x(i-1,5) + w_opt(i)*dt; 
        end
        
        % Objective function to be minimized (Achieved by optimizing dt)
        J = 1*N*dt;
    end
    
end

%% Constraint  Set (State constraints, Control constraints, Path constraints)

function [g,h] = const(dec_vars)

    N = 50;
    dt = dec_vars(end,1);
    a_opt = dec_vars(1:N,1); 
    w_opt = dec_vars(1:N,2);
      x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state1
%     x_init = [10.7 1.95 pi/15 -1 -pi/8];  %%% Initial state2
%     x_init = [10.7 1.95 0 -1 -pi/8];  %%% Initial state3
%     x_init = [10.7 1.95 -pi/9 -1 -pi/8];  %%% Initial state4
    x = zeros(N,5);
    K_dot =zeros(N,1);
    
    % All values are in meters (m)

    Cl = 3.5;                   %%% width of the road
    Lsw = 2;                    %%% width of the parking space
    Lsl = 7;                    %%% Length of the parking space

    L = 2.5;                    %%% Wheel base
    m = 0.7;                    %%% Rear over hang length
    n =0.8;                     %%% Front over hang length

    b = 0.75;                   %%% Half of trackwidth  T = 2b =1.7 m

  % Preallocating corner points vector 
    Ax = zeros(N,1);
    Bx = zeros(N,1);
    Cx = zeros(N,1);
    Dx = zeros(N,1);
    Ay = zeros(N,1);
    By = zeros(N,1);
    Cy = zeros(N,1);
    Dy = zeros(N,1);
        
% Updating the states at each time step using the update equations derived
% by transcribing the vehicle dynamics equations in Continuous time to discrete time
% using Forward Euler approximation
    for i = 1:N

        if i==1
            x(i,1) = x_init(1) + x_init(4)*cos(x_init(3))*dt;
            x(i,2) = x_init(2) + x_init(4)*sin(x_init(3))*dt;
            x(i,3) = x_init(3) + ((x_init(4)*tan(x_init(5)))/L)*dt;
            x(i,4) = x_init(4) + a_opt(i)*dt;
            x(i,5) = x_init(5) + w_opt(i)*dt;
            K_dot(i) = w_opt(i)/(L*cos(x_init(5))^2); % Rate of curvature 
        else
            x(i,1) = x(i-1,1) + x(i-1,4)*cos(x(i-1,3))*dt;
            x(i,2) = x(i-1,2) + x(i-1,4)*sin(x(i-1,3))*dt;
            x(i,3) = x(i-1,3) + ((x(i-1,4)*tan(x(i-1,5)))/L)*dt;
            x(i,4) = x(i-1,4) + a_opt(i)*dt;
            x(i,5) = x(i-1,5) + w_opt(i)*dt; 
            K_dot(i) = w_opt(i)/(L*cos(x(i-1,5))^2); % Rate of curvature 
        end
        
         %Finding the corner points of the vehicle from the states using the
         %following Kinematics equations
         
        Ax(i) = x(i,1) + (L+n)*cos(x(i,3))-b*sin(x(i,3));
        Ay(i) = x(i,2) + (L+n)*sin(x(i,3))+b*cos(x(i,3));

        Bx(i) = x(i,1) + (L+n)*cos(x(i,3))+b*sin(x(i,3));
        By(i) = x(i,2) + (L+n)*sin(x(i,3))-b*cos(x(i,3));

        Cx(i) = x(i,1) - m*cos(x(i,3))+b*sin(x(i,3));
        Cy(i) = x(i,2) - m*sin(x(i,3))-b*cos(x(i,3));

        Dx(i) = x(i,1) - m*cos(x(i,3)) - b*sin(x(i,3));
        Dy(i) = x(i,2) - m*sin(x(i,3)) + b*cos(x(i,3));  
          
    end  

% Pre allocating 'g'

    g = zeros(19*N+1,1);
    g(1) = -N*dt;   %%%% To make the optimization problem well posed
    
%%% Inequality constraints

    for i =1:N
        
        %%% State constraints
        
        g(i+1) = x(i,4)-2;                  %%% Upper bound on velocity
        g(N+i+1) = -x(i,4)-2;               %%% Lower bound on velocity
        g(2*N+i+1) = x(i,5)-pi/6;           %%% Upper bound on steering angle
        g(3*N+i+1) = -x(i,5)-pi/6;          %%% Lower bound on steering angle
        g(4*N+i+1) = x(i,1) -15;            %%% Upper bound on x
        g(5*N+i+1) = -x(i,1)-2;             %%% Lower bound on x
        g(6*N+i+1) = x(i,2)/Cl -1;          %%% Upper bound on y
        g(7*N+i+1) = -x(i,2)/2-1;           %%% Upper bound on y
        
        %%% Control constraints
        
        g(8*N+i+1) = (a_opt(i)/0.25)-1;     %%% Upper bound on optimum acceleration input
        g(9*N+i+1) = (-a_opt(i)/0.75)-1;    %%% Lower bound on optimum acceleration input
        g(21*N+i+1) = w_opt(i)-pi/3;        %%% Upper bound on optimum Steering rate input
        g(22*N+i+1) = -w_opt(i)-pi/3;       %%% Lower bound on optimum Steering rate input
        %%% Path constraints
        
        % To be commented if trying parking 4 scenario
        
%         g(10*N+i+1) = fslot(Ax(i))-Ay(i);   
%         g(11*N+i+1) = Ay(i)-Cl+0.5;
%         g(10*N+i+1) = fslot(Bx(i))-By(i);
%         g(11*N+i+1) = By(i)-Cl+0.5;
%         g(12*N+i+1) = fslot(Cx(i))-Cy(i);
%         g(13*N+i+1) = Cy(i)-Cl+0.5;
%         g(14*N+i+1) = fslot(Dx(i))-Dy(i);
%         g(15*N+i+1) = Dy(i)-Cl+0.5;
%    
        %%% Constraints for smooth trajectory profile
        
        g(16*N+i+1) = K_dot(i)-0.6;         %%% Upper bound on rate of curvature 
        g(17*N+i+1) = -K_dot(i)-0.6;        %%% Lower bound of rate of curvature 

    end

        %%% Constraints for smooth Acceleration profile
    for i = 1:N/2
       g(20*N+i+1) =(a_opt(N/2+i)/0.5)-1;   
    end

%%% Equality constraints (Terminal constraints that make sure the vehicle gets into the desired parking space)
    
    h(1) = x(N,4);          %%% Final velocity should be zero
    h(2) = a_opt(N);        %%% Final acceleration should be zero
    % Parallel parking case
    h(3) = x(N,1)-1;        %%% Desired final position (x)
    h(4) = x(N,2)+1;        %%% Desired final position (y)
%     % Straight parking case
%     h(2) = x(N,1)-1;
%     h(3) = x(N,2)-1;
    
    h(5) = x(N,3);          %%% Desired orientation (theta)

end

%%  Used for Constraining the path

function res = fslot(x)
    Cl = 3.5;                       %%% width of the road
    Lsw = 2;                        %%% width of the parking space
    Lsl = 7.3;                      %%% Length of the parking space
    
    res = -(H(x)+H(x-Lsl))*Lsw;
    
end

%% Heavy side step function (Unit jump function)-- for path constraints

function h =  H(x)

    if (x<0)
        h = 0;
    else
        h = 1;
    end

end

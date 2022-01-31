%% MAGNETIC SATELLITE DYNAMICS %%
%  Contributors: Adyn Miles
%  Description: This code runs continuous magnetic and gravity gradient dynamics.

initial_values = deg2rad([5; 5; 5; 0; 0; 0]);               % Input initial conditions in degrees.


dt = 1;                                                     % Time step
T = 54000;                                                  % Total simulation time
t = 0:dt:T;                                                 % Temporal grid for plotting

% Numerical Integration
options = odeset('RelTol',1e-8,'AbsTol',1e-8);              % Adjust ode45 tolerance
[t, theta] = ode45( @sys, t, initial_values, options);


% Plot Angular Positions
figure;
subplot(3, 1, 1);
plot(t/5400, theta(:,1));
title("Roll Behaviour with Magnetic Dynamics", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 2);
plot(t/5400 , theta(:,2));
title("Pitch Behaviour with Magnetic Dynamics", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 3);
plot(t/5400, theta(:,3));
title("Yaw Behaviour with Magnetic Dynamics", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);


%% State Space Equation with Magnetic
function thetadot = sys(t, theta)
    
    G = 6.672e-11;                                              % Universal gravitational constant
    M_e = 5.98e24;                                              % Mass of Earth
    mu_e = G * M_e;                                             % Gravitational constant of Earth

    I1 = 22;                                                    % Inertia values                              
    I2 = 30;
    I3 = 15;
    
    
    inc = deg2rad(86);                                          % Inclination angle (degrees)

    R_0 = 6871000;                                              % Orbital radius (assume 500km altitude)
    omega_o = sqrt(mu_e/R_0^3);                                 % Natural frequency
    k_1 = (I2 - I3)/I1;                   
    k_3 = (I2 - I1)/I3;
    
    m = [0.1; 0.1; 0.1];                                        % Magnetic dipole moment
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                                              % Magnetic constant 
    
    
    C_io = [1, 0, 0;
            0, cosd(inc), -sind(inc);
            0, sind(inc), cosd(inc)];                           % Rotation from orbital to ECI
        
    C_oi = [cos(omega_o*t), -sin(omega_o*t)*cosd(inc), sin(omega_o*t)
            sin(omega_o*t), cos(omega_o*t)*cosd(inc), -cos(omega_o*t)*sind(inc)
            0, sind(inc), cosd(inc)];                           % Rotation from ECI to perifocal
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];                                            % Rotation to LVLH
    C_bo = zeros(3);
    C_bo(1, 1) = cos(theta(1))*cos(theta(2));
    C_bo(1, 2) = cos(theta(1))*sin(theta(2))*sin(theta(3)) - sin(theta(1))*cos(theta(3));
    C_bo(1, 3) = cos(theta(1))*sin(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(3));
    C_bo(2, 1) = sin(theta(1))*cos(theta(2));
    C_bo(2, 2) = sin(theta(1))*sin(theta(2))*sin(theta(3)) + cos(theta(1))*cos(theta(3));
    C_bo(2, 3) = sin(theta(1))*sin(theta(2))*cos(theta(3)) - cos(theta(1))*sin(theta(3));
    C_bo(3, 1) = -sin(theta(2));
    C_bo(3, 2) = cos(theta(2))*sin(theta(3));
    C_bo(3, 3) = cos(theta(2))*cos(theta(3));                   % Rotation to body-fixed
        
        
    R = C_io * [R_0*cos(omega_o*t); R_0*sin(omega_o*t); 0];     % Orbital position
    X = R(1);
    Y = R(2);
    Z = R(3);
    
    B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];    % Inertial magnetic field
    B_o = C_bo*C_oo*C_oi*B_g;                                   % Body frame magnetic field
    B_ox = [0, B_o(3), -B_o(2);                                 % Skew-symmetric cross matrix
            -B_o(3), 0, B_o(1);
            B_o(2), -B_o(1), 0];

    Ac = zeros(6);                                              % Dynamics matrix
    Ac(1, 4) = 1;
    Ac(2, 5) = 1;
    Ac(3, 6) = 1;
    Ac(4, 1) = -4*k_1*omega_o^2;
    Ac(4, 6) = ((1 - k_1)*omega_o);
    Ac(5, 2) = -(3*(omega_o^2)*(I1 - I3));
    Ac(6, 3) = -(k_3*omega_o^2); 
    Ac(6, 4) = -((1 - k_3)*omega_o);


    Bc_lower = -I_inv*B_ox;

    Bc =[0, 0, 0
        0, 0, 0
        0, 0, 0
        Bc_lower(1,1), Bc_lower(1,2), Bc_lower(1,3);
        Bc_lower(2,1), Bc_lower(2,2), Bc_lower(2,3);
        Bc_lower(3,1), Bc_lower(3,2), Bc_lower(3,3);];          % Build time-varying magnetics matrix

    thetadot = Ac*theta + Bc*m;                                 % Open loop
    fprintf("t = %.1f\n", t);
end
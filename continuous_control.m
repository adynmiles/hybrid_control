%% MAGNETIC SATELLITE CONTROL (PROPORTIONAL-DERIVATIVE) %%
%  Contributors: Adyn Miles
%  Description: This code runs continuous magnetic optimal control using 
%               a proportional-derivative controller.

initial_values = deg2rad([5; 5; 5; 1; 1; 1]);               % Input initial conditions in degrees.



dt = 1;                                                     % Timestep
T = 10800;                                                  % Total time
t = 0:dt:T;                                                 % Time grid for plotting                   

omega = 100;
n = 400;                                                    % Number of turns
d = 0.100;                                                  % Diameter of magnetorquer (m)
A = (pi*d^2)/4;                                             % Area of magnetorquers (m^2)  

options = odeset('RelTol',1e-8,'AbsTol',1e-8);              % Error Tolerance
[t, theta] = ode45( @closed_loop , t, initial_values, options);
m = solve_m(t, theta);                                      % Calculate m(t) for each t solved by the integrator
[Bc, b_o] = solve_Bc(t, theta);

I = m / (n*A);                                              % Calculate input current  
E = (omega)/((n^2) * (A^2))* (T*((m(:,1)'*m(:,1)) + (m(:,2)'*m(:,2)) + (m(:,3)'*m(:,3))));
fprintf("Total Energy Consumed: %dJ\n", E);

torque_mag = zeros(length(t), 3);
for i = 1:length(t)
    m_cross_i = [0, -m(i,3), m(i,2);
                 m(i,3), 0, -m(i,1);
                 -m(i,2), m(i,1), 0];
    torque_mag(i, :) = m_cross_i*b_o(i, :)';
end
    
 
total_torque = sqrt(54000*((torque_mag(:,1)'*torque_mag(:,1) + torque_mag(:,2)'*torque_mag(:,2) + torque_mag(:,3)'*torque_mag(:,3)))/(54000));

% Plot Angular Position
figure;
subplot(3, 1, 1);
plot(t/5400, theta(:,1), 'LineWidth', 2);
title("Roll Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 2);
plot(t/5400, theta(:,2), 'LineWidth', 2);
title("Pitch Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

subplot(3, 1, 3);
plot(t/5400, theta(:,3), 'LineWidth', 2);
title("Yaw Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angle (rad)", 'FontSize', 14);

% Plot Angular Velocity
figure;
subplot(3, 1, 1);
plot(t/5400, theta(:,4), 'LineWidth', 2);
title("Roll Angular Velocity Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

subplot(3, 1, 2);
plot(t/5400, theta(:,5), 'LineWidth', 2);
title("Pitch Angular Velocity Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

subplot(3, 1, 3);
plot(t/5400, theta(:,6), 'LineWidth', 2);
title("Yaw Angular Velocity Behaviour with Magnetic LQR Control", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Angular Velocity (rad/s)", 'FontSize', 14);

% Plot Control Input
figure;
subplot(3, 1, 1); 
plot(t/5400, I(:, 1), 'LineWidth', 2);
title("Control Input Current Required (Roll)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);

subplot(3, 1, 2); 
plot(t/5400, I(:, 2), 'LineWidth', 2);
title("Control Input Current Required (Pitch)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);

subplot(3, 1, 3); 
plot(t/5400, I(:, 3), 'LineWidth', 2);
title("Control Input Current Required (Yaw)", 'FontSize', 14);
xlabel("Number of Orbits", 'FontSize', 14);
ylabel("Current Input (A)", 'FontSize', 14);


%% State Space Equation with Magnetic
function [m] = solve_m(t, theta)
    G = 6.674e-11;                                          % Universal gravitational constant
    M_e = 5.97e24;                                          % Mass of Earth
    mu_e = G * M_e;                                         % Gravitational constant of Earth

    I1 = 15;                                                % Inertia values                            
    I2 = 30;
    I3 = 22;
    
    %K_p = 200000;                                           % Proportional gain constant
    %K_d = 100000000;                                        % Derivative gain constant
    K_p = 1000000;
    K_d = 500000000;
    
    inc = deg2rad(86);                                      % Inclination Angle (degrees) 

    R_0 = 6871000;                                          % Orbital radius (assuming 500km altitude)
    omega_o = sqrt(mu_e/R_0^3);                             % Natural frequency
    k_1 = (I2 - I3)/I1; 
    k_3 = (I2 - I1)/I3;
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                                          % Magnetic constant 
    
    
    C_io = [1, 0, 0;
            0, cos(inc), -sin(inc);
            0, sin(inc), cos(inc)];                         % Rotation from orbital to ECI
        
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];                                        % Rotation to LVLH
        
    m = zeros(length(t), 3);
    C_bo = zeros(3);
        
    
    
    for i = 1:length(t)
        C_oi = [cos(omega_o*t(i)), -sin(omega_o*t(i))*cos(inc), sin(omega_o*t(i));
                sin(omega_o*t(i)), cos(omega_o*t(i))*cos(inc), -cos(omega_o*t(i))*sin(inc);
                0, sin(inc), cos(inc)];                     % Rotation from ECI to perifocal
        
        C_bo(1, 1) = cos(theta(i, 1))*cos(theta(i, 2));
        C_bo(1, 2) = cos(theta(i, 1))*sin(theta(i, 2))*sin(theta(i, 3)) - sin(theta(i, 1))*cos(theta(i, 3));
        C_bo(1, 3) = cos(theta(i, 1))*sin(theta(i, 2))*cos(theta(i, 3)) + sin(theta(i, 1))*sin(theta(i, 3));
        C_bo(2, 1) = sin(theta(i, 1))*cos(theta(i, 2));
        C_bo(2, 2) = sin(theta(i, 1))*sin(theta(i, 2))*sin(theta(i, 3)) + cos(theta(i, 1))*cos(theta(i, 3));
        C_bo(2, 3) = sin(theta(i, 1))*sin(theta(i, 2))*cos(theta(i, 3)) - cos(theta(i, 1))*sin(theta(i, 3));
        C_bo(3, 1) = -sin(theta(i, 2));
        C_bo(3, 2) = cos(theta(i, 2))*sin(theta(i, 3));
        C_bo(3, 3) = cos(theta(i, 2))*cos(theta(i, 3));     % Rotation to body-fixed  

        % Orbital position
        R = C_io * [R_0*cos(omega_o*t(i)); R_0*sin(omega_o*t(i)); 0];
        X = R(1);
        Y = R(2);
        Z = R(3);
        
        % Inertial magnetic field
        B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];
        B_o = C_bo*C_oo*C_oi*B_g;                           % Body-fixed magnetic field
        B_ox = [0, -B_o(3), B_o(2);                         % Body-fixed cross matrix
                B_o(3), 0, -B_o(1);
                -B_o(2), B_o(1), 0];
            
        K = [K_p, 0, 0, K_d, 0, 0;                          % Constant gain matrix
            0, K_p, 0, 0, K_d, 0;
            0, 0, K_p, 0, 0, K_d];
        Kc = -B_ox * K;                                     % Time-varying gain matrix
        
        m(i, :) = Kc*theta(i, :)';                          % Calculation of magnetic dipole moment
        
    end
end

function [thetadot] = closed_loop(t, theta)
    G = 6.674e-11;                      
    M_e = 5.97e24;                      
    mu_e = G * M_e;                     

    I1 = 15;                            
    I2 = 30;
    I3 = 22;
    
    %K_p = 200000; 
    %K_d = 100000000;
    K_p = 1000000;
    K_d = 500000000;
    
    inc = deg2rad(86);                             

    R_0 = 6871000;                      
    omega_o = sqrt(mu_e/R_0^3);
    k_1 = (I2 - I3)/I1; 
    k_3 = (I2 - I1)/I3;
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                       
    
    
    C_io = [1, 0, 0;
            0, cos(inc), -sin(inc);
            0, sin(inc), cos(inc)];
        
    C_oi = [cos(omega_o*t), -sin(omega_o*t)*cos(inc), sin(omega_o*t);
            sin(omega_o*t), cos(omega_o*t)*cos(inc), -cos(omega_o*t)*sin(inc);
            0, sin(inc), cos(inc)];
        
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];
    C_bo = zeros(3);
    C_bo(1, 1) = cos(theta(1))*cos(theta(2));
    C_bo(1, 2) = cos(theta(1))*sin(theta(2))*sin(theta(3)) - sin(theta(1))*cos(theta(3));
    C_bo(1, 3) = cos(theta(1))*sin(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(3));
    C_bo(2, 1) = sin(theta(1))*cos(theta(2));
    C_bo(2, 2) = sin(theta(1))*sin(theta(2))*sin(theta(3)) + cos(theta(1))*cos(theta(3));
    C_bo(2, 3) = sin(theta(1))*sin(theta(2))*cos(theta(3)) - cos(theta(1))*sin(theta(3));
    C_bo(3, 1) = -sin(theta(2));
    C_bo(3, 2) = cos(theta(2))*sin(theta(3));
    C_bo(3, 3) = cos(theta(2))*cos(theta(3));
        
        
    R = C_io * [R_0*cos(omega_o*t); R_0*sin(omega_o*t); 0];
    X = R(1);
    Y = R(2);
    Z = R(3);
    
    B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];
    B_o = C_bo*C_oo*C_oi*B_g;
    B_ox = [0, -B_o(3), B_o(2);
            B_o(3), 0, -B_o(1);
            -B_o(2), B_o(1), 0];

    Ac = zeros(6);
    Ac(1, 4) = 1;
    Ac(2, 5) = 1;
    Ac(3, 6) = 1;
    Ac(4, 1) = -4*k_1*omega_o^2;
    Ac(4, 6) = ((1 - k_1)*omega_o);
    Ac(5, 2) = -(3*(omega_o^2)*(I1 - I3));
    Ac(6, 3) = -(k_3*omega_o^2); 
    Ac(6, 4) = -((1 - k_3)*omega_o);


    Bc_lower = -I_inv*B_ox;

    Bc =[0, 0, 0;
        0, 0, 0;
        0, 0, 0;
        Bc_lower(1,1), Bc_lower(1,2), Bc_lower(1,3);
        Bc_lower(2,1), Bc_lower(2,2), Bc_lower(2,3);
        Bc_lower(3,1), Bc_lower(3,2), Bc_lower(3,3);];
    
    K = [K_p, 0, 0, K_d, 0, 0;
        0, K_p, 0, 0, K_d, 0;
        0, 0, K_p, 0, 0, K_d];
    Kc = -B_ox * K; 
    thetadot = (Ac - (Bc*Kc))*theta;                        % Closed-loop feedback
    
    fprintf("t = %.1f\n", t);
end

function [Bc, B_o] = solve_Bc(t, theta)
    G = 6.674e-11;                      
    M_e = 5.97e24;                      
    mu_e = G * M_e;                     

    I1 = 15;                            
    I2 = 30;
    I3 = 22;
    
    K_p = 1000000;
    K_d = 500000000;
    
    inc = deg2rad(86);                             

    R_0 = 6871000;                      
    omega_o = sqrt(mu_e/R_0^3);
    k_1 = (I2 - I3)/I1; 
    k_3 = (I2 - I1)/I3;
    
    I_inv = [-1/I1, 0, 0;
            0, -1/I2, 0;
            0, 0, -1/I3];
        
    B_0 = -8.0e15;                      
    
    
    C_io = [1, 0, 0;
            0, cos(inc), -sin(inc);
            0, sin(inc), cos(inc)];
    C_oo = [0 1 0; 
            0 0 -1;
            -1 0 0];
        
    Bc = zeros(length(t), 6, 3);
    B_o = zeros(length(t), 3);
    
    for i = 1:length(t)
        C_oi = [cos(omega_o*t(i)), -sin(omega_o*t(i))*cos(inc), sin(omega_o*t(i));
                sin(omega_o*t(i)), cos(omega_o*t(i))*cos(inc), -cos(omega_o*t(i))*sin(inc);
                0, sin(inc), cos(inc)];
        
        C_bo = eye(3);
        
        R = C_io * [R_0*cos(omega_o*t(i)); R_0*sin(omega_o*t(i)); 0];
        
        X = R(1);
        Y = R(2);
        Z = R(3);

        B_g = (B_0 / R_0^5) * [3*X*Z; 3*Y*Z; 2*Z^2 - X^2 - Y^2];
        B_o(i, :) = C_bo*C_oo*C_oi*B_g;
        B_ox = [0, -B_o(i, 3), B_o(i, 2);
                B_o(i, 3), 0, -B_o(i, 1);
                -B_o(i, 2), B_o(i, 1), 0];
        Bc_lower = -I_inv*B_ox;

        Bc(i, :, :) =[0, 0, 0;
                        0, 0, 0;
                        0, 0, 0;
                        Bc_lower(1,1), Bc_lower(1,2), Bc_lower(1,3);
                        Bc_lower(2,1), Bc_lower(2,2), Bc_lower(2,3);
                        Bc_lower(3,1), Bc_lower(3,2), Bc_lower(3,3);];
    end
end

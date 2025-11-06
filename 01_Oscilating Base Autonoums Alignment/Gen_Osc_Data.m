clc;
clearvars;
close all;
rng(5);   %% for regeneration
%%
dt  = 0.02;
d2r = pi/180;
r2d = 1/d2r;
Lat = 20.0*d2r;
H   = 0;
ge  = 9.78049;
Wie = 7.2921158e-5;
Re  = 6378137;
beta  = 0.005317;
alpha = 0.000007;
g = ge*Re^2/(Re+H)^2*(1 + beta*sin(Lat)^2 + alpha*sin(2*Lat)^2);
gen_data = 1;
if gen_data == 1
    t_data = 15*60; % sec
    dt_ref = 0.02;
    t_ref  = (0:dt_ref:t_data)';

    osc_flag = 1; % Enable/Disable Angular Oscillations
    acc_flag = 1; % Enable/Disable Linear Oscillations

    A_r = 10*osc_flag;
    A_p = 7*osc_flag;
    A_h = 5*osc_flag; % Amplitude
    T_r =          6;
    T_p =          5;
    T_h =          7; % Time Period
    P_r =  pi/7;
    P_p =  pi/4;
    P_h =  pi/3; % Phase
    P_r =  pi/7; %2*pi*rand;
    P_p =  pi/4; %2*pi*rand;
    P_h =  pi/3; %2*pi*rand; % Phase

    roll_ref    = (5 + A_r*cos(2*pi*t_ref/T_r + P_r))*d2r;
    pitch_ref   = (10 + A_p*cos(2*pi*t_ref/T_p + P_p))*d2r;
    heading_ref = (30 + A_h*cos(2*pi*t_ref/T_h + P_h))*d2r;

    roll_dot  = -A_r*2*pi/T_r*sin(2*pi*t_ref/T_r + P_r)*d2r;
    pitch_dot = -A_p*2*pi/T_p*sin(2*pi*t_ref/T_p + P_p)*d2r;
    head_dot  = -A_h*2*pi/T_h*sin(2*pi*t_ref/T_h + P_h)*d2r;

    P = pitch_dot.*cos(roll_ref) - head_dot.*sin(roll_ref).*cos(pitch_ref);
    Q = roll_dot + head_dot.*sin(pitch_ref);
    R = pitch_dot.*sin(roll_ref) + head_dot.*cos(roll_ref).*cos(pitch_ref);

    A_x = 0.3*acc_flag;
    A_y = 0.03*acc_flag;
    A_z = 0.02*acc_flag; % Amplitude
    T_x =          8;
    T_y =          6;
    T_z =          7; % Time Period
    % P_x =  pi;
    % P_y =  pi;
    % P_z =  pi; % Phase
    P_x =  2*pi*rand;
    P_y =  2*pi*rand;
    P_z =  2*pi*rand; % Phase

    X  = A_x*cos(2*pi*t_ref/T_x + P_x);
    Y  = A_y*cos(2*pi*t_ref/T_y + P_y);
    Z  = A_z*cos(2*pi*t_ref/T_z + P_z);

    Vx = -A_x*(2*pi/T_x)*sin(2*pi*t_ref/T_x + P_x);
    Vy = -A_y*(2*pi/T_y)*sin(2*pi*t_ref/T_y + P_y);
    Vz = -A_z*(2*pi/T_z)*sin(2*pi*t_ref/T_z + P_z);

    Ax = -A_x*(2*pi/T_x)^2*cos(2*pi*t_ref/T_x + P_x);
    Ay = -A_y*(2*pi/T_y)^2*cos(2*pi*t_ref/T_y + P_y);
    Az = -A_z*(2*pi/T_z)^2*cos(2*pi*t_ref/T_z + P_z);

    save reference_Osc t_ref heading_ref pitch_ref roll_ref Ax Ay Az Vx Vy Vz X Y Z P Q R dt_ref r2d dt;

    mg2mpsps = 9.8/1000;
    dph2rps = pi/180/3600;
    b_flag = 1; % Enable/Disable Sensor Bias
    n_flag = 0; % Enable/Disable Sensor Noise
    acc_bias   = 0.1*mg2mpsps*b_flag;
    gyro_bias  = 0.01*dph2rps*b_flag;
    acc_noise  = 0.05*mg2mpsps*n_flag;
    gyro_noise = 0.05*dph2rps*n_flag;

    %% Write IMU Data
    fid = fopen('IMU_Data_Osc.dat','w');
    %    for i = 1 : dt/dt_ref:  length(t_ref)-dt/dt_ref
    for i = 1 : length(t_ref)


        Cn2b = angle2dcm(heading_ref(i),pitch_ref(i),roll_ref(i),'zxy');
        Ab  = Cn2b*[0;0;g] + [Ax(i);Ay(i);Az(i)];
        Wb  = [P(i);Q(i);R(i)] + Cn2b*[0;Wie*cos(Lat);Wie*sin(Lat)];

        dvx = Ab(1)*dt + acc_bias*dt + acc_noise*randn*dt;
        dvy = Ab(2)*dt + acc_bias*dt + acc_noise*randn*dt;
        dvz = Ab(3)*dt + acc_bias*dt + acc_noise*randn*dt;

        dax = Wb(1)*dt + gyro_bias*dt + gyro_noise*randn*dt;
        day = Wb(2)*dt + gyro_bias*dt + gyro_noise*randn*dt;
        daz = Wb(3)*dt + gyro_bias*dt + gyro_noise*randn*dt;

        fprintf(fid,'%10.5f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n',t_ref(i),dax,day,daz,dvx,dvy,dvz);
    end
    fclose('all');
end

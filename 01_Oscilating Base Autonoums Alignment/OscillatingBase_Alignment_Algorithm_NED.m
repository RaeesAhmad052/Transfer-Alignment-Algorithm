clc;clear all;close all;

dt  = 0.02;
d2r = pi/180;
r2d = 1/d2r;
Lat = 33.0*d2r;
H   = 0;
ge  = 9.78049;
Wie = 7.2921158e-5;
Re  = 6378137;

beta  = 0.005317;
alpha  = 0.000007;
g = ge*Re*Re/(Re+H)/(Re+H)*(1.+ beta*sin(Lat)*sin(Lat) + alpha*sin(2*Lat)*sin(2*Lat));

%% Read IMU Date Populate in body frame x-forward, y-rightward, z-downward
IMU = load ('IMU_Data_Osc.dat');
Time= IMU(:,1);
Abx = IMU(:,6)/dt;
Aby = IMU(:,5)/dt;
Abz =-IMU(:,7)/dt;
Wbx = IMU(:,3)/dt;
Wby = IMU(:,2)/dt;
Wbz =-IMU(:,4)/dt;

%% Alignment Algorithm Initilization
fid = fopen('OBA_sim_ca.dat','w');

Vi  = [0;0;0];
Vib = [0;0;0];

Roll  = 0;Pitch  = 0;Heading  = 0;
R  = 0;P  = 0;H  = 0;

Cib2i  = zeros(3,3);

%% Initialize the Transformation from  inerital Body frame to Current Body frame
qib2b0t1 = 1;
qib2b1t1 = 0;
qib2b2t1 = 0;
qib2b3t1 = 0;

%% Transformation from ECEF frame to ENU frame
Ce2n  = [-sin(Lat) 0 cos(Lat);0 1 0;-cos(Lat) 0 -sin(Lat)];

%% Alignment Loop
%aligntime = min([Time(end) 600]);
aligntime = 600;
t1 = 450;
t2 = 550;
t = 0;
for i = 1 : length(IMU)
    t  = t + dt;
    
    Ab = [Abx(i);Aby(i);Abz(i)];
    Wb = [Wbx(i);Wby(i);Wbz(i)];

    %% Transformation from Inertial Body frame to current Body frame
    phix = Wb(1)*dt;
    phiy = Wb(2)*dt;
    phiz = Wb(3)*dt;
    phi  = phix*phix + phiy*phiy + phiz*phiz;
    
    % f3  = sin(sqrt(phi)/2)/sqrt(phi);  
    % f4  = cos(sqrt(phi)/2);            

    f3  = 1/2 - phi/48 + phi*phi/3840;
    f4  = 1 - phi/8 + phi*phi/384;
    
    qib2b0 = f4*qib2b0t1 - f3*( phix*qib2b1t1 + phiy*qib2b2t1 + phiz*qib2b3t1);
    qib2b1 = f4*qib2b1t1 + f3*( phix*qib2b0t1 - phiy*qib2b3t1 + phiz*qib2b2t1);
    qib2b2 = f4*qib2b2t1 + f3*(-phiz*qib2b1t1 + phix*qib2b3t1 + phiy*qib2b0t1);
    qib2b3 = f4*qib2b3t1 + f3*( phiz*qib2b0t1 - phix*qib2b2t1 + phiy*qib2b1t1);
    
    dqib2b = 1 - (qib2b0*qib2b0 + qib2b1*qib2b1 + qib2b2*qib2b2 + qib2b3*qib2b3);
    dqib2b = 1 + dqib2b/2;
    
    qib2b0 = qib2b0*dqib2b;
    qib2b1 = qib2b1*dqib2b;
    qib2b2 = qib2b2*dqib2b;
    qib2b3 = qib2b3*dqib2b;
    
    qib2b0t1 = qib2b0;
    qib2b1t1 = qib2b1;
    qib2b2t1 = qib2b2;
    qib2b3t1 = qib2b3;
    
    Cib2b(1,1) = 1 - 2*(qib2b2t1*qib2b2t1 + qib2b3t1*qib2b3t1);
    Cib2b(1,2) =     2*(qib2b1t1*qib2b2t1 + qib2b0t1*qib2b3t1);
    Cib2b(1,3) =     2*(qib2b1t1*qib2b3t1 - qib2b0t1*qib2b2t1);
    Cib2b(2,1) =     2*(qib2b1t1*qib2b2t1 - qib2b0t1*qib2b3t1);
    Cib2b(2,2) = 1 - 2*(qib2b1t1*qib2b1t1 + qib2b3t1*qib2b3t1);
    Cib2b(2,3) =     2*(qib2b2t1*qib2b3t1 + qib2b0t1*qib2b1t1);
    Cib2b(3,1) =     2*(qib2b1t1*qib2b3t1 + qib2b0t1*qib2b2t1);
    Cib2b(3,2) =     2*(qib2b2t1*qib2b3t1 - qib2b0t1*qib2b1t1);
    Cib2b(3,3) = 1 - 2*(qib2b1t1*qib2b1t1 + qib2b2t1*qib2b2t1);
    

    Cb2ib = transpose(Cib2b);

    %% Transformation from ECI frame to ECEF frame
    Ci2e = [cos(Wie*t) sin(Wie*t) 0;-sin(Wie*t) cos(Wie*t) 0;0 0 1];
    
    %% Velocity Computation in ECI and Inertial Body Frame
    Vi   = [g*cos(Lat)*sin(Wie*t)/Wie;g*cos(Lat)*(1-cos(Wie*t))/Wie;t*g*sin(Lat)];
    dVib = Cib2b'*Ab*dt;
    dWib = Cib2b'*Wb*dt;
    Vib  = Vib + dVib + 0.5*cross(dWib,dVib);
    
    %% Transformation from Inertial body to ECI frame
    if round(t*100)/100 == t1
        Vi1  = Vi;
        Vib1 = Vib;
    elseif t >= t2
        Vi2  = Vi;
        Vib2 = Vib;
        V1 = [Vi1' ;Vi2' ;cross(Vi1 ,Vi2)' ];
        V2 = [Vib1';Vib2';cross(Vib1,Vib2)'];
        Cib2i = V1\V2;
        
        %% Transformation from body frame to Navigation (ENU) 
        Cb2n = Ce2n*Ci2e*Cib2i*Cb2ib;
                

        %% Attitude Computation      
        
        Heading  = atan2(Cb2n(2,1),Cb2n(1,1));
        Roll  = atan2(Cb2n(3,2),Cb2n(3,3));
        Pitch = atan2( -Cb2n(3,1),sqrt(Cb2n(3,2)^2+Cb2n(3,3)^2));

    end
    
    %% Save Parameters
    fprintf(fid,'%07.03f\t',Time(i));
    
    fprintf(fid,'%12.12f\t',Roll*r2d);
    fprintf(fid,'%12.12f\t',Pitch*r2d);
    fprintf(fid,'%12.12f\t',-Heading*r2d);


    fprintf(fid,'\n');
end
fclose(fid);
load OBA_sim_ca.dat;
%load OBA_sim_ca;
load reference_Osc;
%% PLOTS
close all;
h=figure(1);set(h,'name','Attitude');
ax(1)=subplot(311);plot(t_ref,   roll_ref*r2d,'k.-',OBA_sim_ca(:,1),OBA_sim_ca(:,2),'r.-');ylabel('Roll');title('Attitude (deg)');grid;legend('Reference','OBA Align Algo','Location','best');
ax(2)=subplot(312);plot(t_ref,  pitch_ref*r2d,'k.-',OBA_sim_ca(:,1),OBA_sim_ca(:,3),'r.-');ylabel('Pitch');grid;
ax(3)=subplot(313);plot(t_ref,heading_ref*r2d,'k.-',OBA_sim_ca(:,1),OBA_sim_ca(:,4),'r.-');ylabel('Heading');xlabel('Alignment Time (sec)');grid;shg;
linkaxes(ax,'x');grid on;shg;
%%
index_OBA = t2/dt;
index_ref = t2/dt_ref;
h=figure(2);set(h,'name','Attitude Error');
subplot(311);plot(t_ref(index_ref:dt/dt_ref:end),  (-roll_ref(index_ref:dt/dt_ref:end)*r2d    +   OBA_sim_ca(index_OBA:end,2)),'r.-');ylabel('Roll');title('Attitude Error (deg)');grid;xlim([t2+1 t]);
subplot(312);plot(t_ref(index_ref:dt/dt_ref:end),  (-pitch_ref(index_ref:dt/dt_ref:end)*r2d   +   OBA_sim_ca(index_OBA:end,3)),'r.-');ylabel('Pitch');grid;xlim([t2+1 t]);
subplot(313);plot(t_ref(index_ref:dt/dt_ref:end),  (-heading_ref(index_ref:dt/dt_ref:end)*r2d + OBA_sim_ca(index_OBA:end,4)),'r.-');ylabel('Heading');xlabel('Alignment Time (sec)');grid;shg;xlim([t2+1 t]);

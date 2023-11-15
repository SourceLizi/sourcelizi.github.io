global G lambda epslion v_T
G=9.7833; %m^2/s
% R=(42.5/2)*1E-3;  %42mm,m
% m = 43*1E-3; %42mm, kg
R=(16.8/2)*1E-3; %17mm,m
m = 3.2*1E-3; %17mm, kg

theta_armor=75/360*(2*pi); %rad
armor_len=125E-3; %m
epslion = 1E-5;

%F=C_d/2*S*ρ*v^2, C_d=0.47
lambda = (0.47/2*1.205*pi*R*R)/m;
v_T=sqrt(G/lambda); %terminal speed

V0=15;
t=[1e-5:0.01:1]; % shorter time, less iteration for solver
theta0 = [0:0.2:90]/360*(2*pi); %rad 
[x,z] = drag2_peak(V0*cos(theta0),V0*sin(theta0));
[x_neg,z_neg] = drag2_trajectory(t,V0,0);
hold on;
plot(x,z)
plot(x_neg,z_neg)
grid on;
ylabel('z/m')
xlabel('x/m')
hold off;

function [x,z] = nodrag_trajectory(t,vx0,vz0)
    global G;
    x = vx0*t;
    z = vz0*t-G*t.*t/2;
end

function [x,z] = drag2_trajectory(t,vx0,vz0)
    global G lambda;
    x=log(1+lambda*vx0*t)/lambda;
    beta = atan(vz0*sqrt(lambda/G));
    if vz0 > 0
        z=log(cos(beta-sqrt(G*lambda)*t)/cos(beta))/lambda.*(t<=beta/sqrt(G*lambda)) ...
        -log(cos(beta)*cosh(sqrt(G*lambda)*t-beta))/lambda.*(t>beta/sqrt(G*lambda));
    else
        z = -log(cosh(sqrt(G*lambda)*t)-tan(beta)*sinh(sqrt(G*lambda)*t))/lambda;
    end
end

function [x,z] = drag2_peak(vx0,vz0)
    global G lambda;
    beta = atan(vz0*sqrt(lambda/G));
    x=log(1+lambda*vx0.*beta/sqrt(G*lambda))/lambda;
    z=-log(cos(atan(vz0*sqrt(lambda/G))))/lambda;
end
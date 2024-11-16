% Modeling for Lift Augmentation and Thrust Modeling
% Author: Gabriel Rodriguez
% VerdeCommute Senior Design
clc; clear; close all;

%% Thrust and Lift Augmentation Model Inputs

%------------------------ Propeller Configuration ------------------------%

N_hlp = 8;                  % number of high lift props (hlp)
Kp_hlp = 0.30;              % hlp scale based correction factor
N_tip = 2;                  % propellers for thrust production
Kp_tip = 0.75;              % tip props scale based correction factor
a_desired = 0.6;            % axial velocity induction factor
u_hlp = 0.75;               % distance upstream of wing leading edge
u_tip = 1.00;               % distance upstream of wing leading edge
d_spinner_hlp = 0.5;        % ft
d_spinner_tip = 0.65;       % ft
num_blades_hlp = 5;         % prop blades
num_blades_tip = 5;         % prop blades

%---------------- Engine, Motor, and Thrust Requirements -----------------%

P_eng = 556;                % bhp
T_req_TO = 2000;            % lbf
motor_RPM_hlp = 4500;       % rpm
motor_RPM_tip = 3500;       % rpm

%------------------------ Aircraft Configuration -------------------------%

b_wing = 56.14;             % ft
fus_width = 6;              % ft
V_TO = 60;                  % KTAS
V_cruise = 150;             % KTAS
c_root = 6.27;              % ft
c_tip = 2.51;               % ft
serv_ceiling = 10001;       % ft
altitude = 0;               % ft

%% Additional Constants from Inputs

P_eng_tot = 550*P_eng;          % Engine Power in ft*lbf/sec
N_tot = N_tip + N_hlp;          % Total number of props
R = 1716;                                               % ftlbf/slugmoleR
gamma = 1.4;                                            % 1
deltaT = 0;                                             % degrees R
hs = linspace(0, serv_ceiling, serv_ceiling);                         % ft
thetas = 1-0.0000068756.*hs;                            % 1
ps = 2116.22.*thetas.^5.2561;                           % lbf/ft^2
Ts_tip = 518.67 * thetas + deltaT;                          % degrees R
rhos = ps ./ (R.*Ts_tip);                                   % slugs/ft^3
rho_des = rhos(1,1);
rho_model = rhos(1,altitude+1);
A_s = pi*d_spinner_hlp^2/4;                             % ft^2
space = (b_wing-fus_width)/N_hlp;                       % ft

%% Using Momentum Theory for  High Lift Prop Sizing

% Distributing 1/8 of engine power to the HLP's
P_hlp = 14.08073*550;                               % ft*lbf/sec per prop
T_hlp_req = (T_req_TO * (1/8))/N_hlp;               % lbf per prop

% Determining Polynomial based on Equation 15-81 from Gudmundsson
A = 1;
B = -(3*A_s + T_hlp_req^3/Kp_hlp^3/P_hlp^2/2/rho_des);
C = 3*A_s^2;
D = -A_s^3;
coeffs = [A B C D];

% Solving for High Lift Prop Radius 
A_p_hlp = real(max(roots(coeffs)));                  % ft^2
r_p_hlp = (sqrt(4*A_p_hlp/pi)/2)/0.85;               % ft, scale factor

%% Using Momentum Theory for Wing Tip Prop Sizing

% Distributing 3/4 of engine power to the HLP's
P_tip = (P_eng_tot-P_hlp*N_hlp)/N_tip;               % ft*lbf/sec per prop
T_tip_req = 600;                                     % lbf per prop

% Determining Polynomial based on Equation 15-81 from Gudmundsson
A = 1;
B = -(3*A_s + T_tip_req^3/Kp_tip^3/P_tip^2/2/rho_des);
C = 3*A_s^2;
D = -A_s^3;
coeffs = [A B C D];

% desired propeller area for 
A_p_tip = real(max(roots(coeffs)));                 % ft^2
w_p_tip = sqrt(T_tip_req/2/rho_des/A_p_tip);      % ft/s

% Solving for Tip Prop Radius and Pitch, acconting for eta expected
r_p_tip = (sqrt(4*A_p_tip/pi)/2)*1.9;               % ft, scale factot
P_G_tip = 1251*((V_TO+V_cruise)/2/motor_RPM_tip);   % Geo Pitch
Beta_tip = atand(P_G_tip/2/pi/r_p_tip);             % Expected Pitch Angle

%% Fixed Pitch Thrust Model - Tip Propellers

eta_opt = 0.75;
A_p_tip = pi*r_p_tip^2;
A_spinner_tip = pi*d_spinner_tip^2/4;
T_static = Kp_tip*(P_tip)^(2/3)* ...
    (2*rho_model*A_p_tip)^(1/3)*(1-A_spinner_tip/A_p_tip);
V_opt = (V_TO+V_cruise)/2.5*1.68781;
T_opt = eta_opt*P_tip/V_opt;
V_opt = (V_TO+V_cruise)/2;
A_Vs = (T_static-2*T_opt)/V_opt^2;
B_Vs = (3*T_opt-2*T_static)/V_opt;
C_Vs = T_static;
Vs = linspace(0, V_cruise + 50);
Ts_tip = A_Vs.*Vs.^2 + B_Vs.*Vs + C_Vs;
T_total_tip = 2*(A_Vs.*60^2 + B_Vs*60 + C_Vs);

%% Designing the High Lift Propellers

% airfoil cls and cds txt file
file_path_hlp = 'xf-mh114-il-500000.txt';
airfoil_data_hlp = readtable(file_path_hlp, 'FileType', 'text',...
    'HeaderLines', 12);

% cls and cds for each angle of attack for HLP (117 data points)
aoas_hlp = airfoil_data_hlp{:, 1};
cls_hlp = airfoil_data_hlp{:, 2};
cds_hlp = airfoil_data_hlp{:, 3};

% Load data from Maxewell X-57 chord distribution digitized plot
chord_dis_data = readtable('x57_HLP_chord_distribution.csv');
r_R = chord_dis_data.x; % location normalized by prop radius
c_R = chord_dis_data.y; % chord normalized by prop radius

% Setting axial induction factor distribution
B = num_blades_hlp;                                 % blades
Omega = motor_RPM_hlp*2*pi/60;                      % rad/s
V_inf = V_TO*1.68781;                               % ft/s
fidelity = 100;                                     % elements
rs_hlp = linspace(d_spinner_hlp/2, r_p_hlp, fidelity);              % ft
helixs = atan(V_inf/Omega./rs_hlp);                                 % rad
f = B/2*(r_p_hlp-rs_hlp)./(rs_hlp.*sin(helixs));
f(length(f)) = f(length(f)-1);
F = 2/pi*acos(exp(-f));                                 % tip loss factor
as = a_desired./F; 
as_prime = (1-sqrt(1-(4*V_inf^2.*(1+as).*as)/Omega^2./rs_hlp.^2))/2; 
as_prime = real(as_prime);
V_Es = sqrt((V_inf*(1+as)).^2+(Omega.*rs_hlp).^2);                  % ft/s

% chord from digitized plot
r_Rs = linspace(min(r_R), max(r_R), fidelity);
c_R_spline = spline(r_R, c_R, r_Rs);
cs_hlp = c_R_spline.*r_p_hlp;

% determining optimal angles of attack for each element
aoa_induced = asin(V_inf*as./sqrt(V_inf^2+(Omega.*rs_hlp).^2)); % rad
dTs = 4*pi.*rs_hlp*rho_des*V_inf^2.*(1+as).*as.*F;            % lbf
dTs_AOAs = zeros(length(aoas_hlp), length(rs_hlp));             % lbf
dTs_diff = zeros(length(aoas_hlp), length(rs_hlp));             % lbf
min_dT_diffs = zeros(1, length(rs_hlp)); 
aoa_inds = zeros(1, length(rs_hlp));
V_ps = zeros(1, length(rs_hlp));
for els=1:length(rs_hlp)
    for aoa=1:length(aoas_hlp)
        dTs_AOAs(aoa,els) = B*cs_hlp(els)*0.5*rho_des.*V_Es(els).^2 ...
            *(cls_hlp(aoa)*cos(helixs(els))-cds_hlp(aoa)*cos(helixs(els)));
        dTs_diff(aoa, els) = abs(dTs(1,els) - dTs_AOAs(aoa, els));
    end
    [min_dT_diffs(els), aoa_inds(els)] = min(dTs_diff(:,els));
end

% obtaining final pitch distribution
aoa_opts_hlp = aoas_hlp(aoa_inds);
cl_dis_hlp = cls_hlp(aoa_inds);
cd_dis_hlp = cds_hlp(aoa_inds);
pitch_rad = helixs.' + aoa_induced.' + aoa_opts_hlp*pi/180;
pitch_dis = pitch_rad*180/pi;

% smoothing out tip pitch distribution
for i = 2:length(pitch_dis)
    slope = (pitch_dis(i)-pitch_dis(i-1))/(rs_hlp(i)-rs_hlp(i-1));
    if slope > 20
        pitch_dis(i) = pitch_dis(i-1);
    end
end

% Fit a smoothing spline to the pitch distribution
smoothing_factor = 0.99; % Adjust between 0 (very smooth) to 1 (exact fit)
spline_fit = fit(rs_hlp(:) / r_p_hlp, pitch_dis(:), 'smoothingspline', 'SmoothingParam', smoothing_factor);
pitch_dis_spline = feval(spline_fit, rs_hlp / r_p_hlp);

% solving for the thrust
dTs_actual = dTs_AOAs(aoa_inds);
T_per_prop_hlp = trapz(rs_hlp, dTs_actual);
T_total_hlp = T_per_prop_hlp*N_hlp;
eta_p_hlp = T_total_hlp*V_inf/550/P_hlp;

% solving for total thrust
T_total = T_total_hlp + T_total_tip;

%% Fixed Pitch Thrust Model - High Lift Propellers

n = motor_RPM_hlp/60;
C_T = T_per_prop_hlp/rho_des/n^2/(r_p_hlp*2)^4;
C_P = P_hlp/rho_des/n^3/(r_p_hlp*2)^5;
J = V_inf/n/(r_p_hlp*2);
eta_hlp = J * C_T/C_P;
A_spinner_hlp = pi*d_spinner_hlp^2/4;
T_static = Kp_hlp*(P_hlp)^(2/3)* ...
     (2*rho_model*A_p_hlp)^(1/3)*(1-A_spinner_hlp/A_p_hlp);
V_opt = V_inf;
T_opt = eta_hlp*P_hlp/V_opt;
V_opt = V_TO;
A_Vs_hlp = (T_static-2*T_opt)/V_opt^2;
B_Vs_hlp = (3*T_opt-2*T_static)/V_opt;
C_Vs_hlp = T_static;
Ts_hlp = A_Vs_hlp.*Vs.^2 + B_Vs_hlp.*Vs + C_Vs_hlp;
pos_indices = Ts_hlp > 0; 
Ts_hlp = Ts_hlp(pos_indices);  
Vs_hlp = Vs(pos_indices); 

%% Total Fixed Pitch Thrust Model

T_tots = Ts_tip*N_tip;
for i = 1:length(Vs_hlp)
    T_tots(1,i) = Ts_hlp(1,i)*N_hlp + Ts_tip(1+i)*N_tip;
end

% Getting polynomial
mean_Vs = mean(Vs);
std_Vs = std(Vs);
Vs_scaled = (Vs - mean_Vs)/std_Vs;
Ts_tots_func = polyfit(Vs_scaled, T_tots, 10);
T_fit = polyval(Ts_tots_func, Vs_scaled);

%% Velocity Profile aft of High Lift Props

for els = 1:length(rs_hlp)
    V_ps(1,els) = v_induced(V_inf, Omega, rs_hlp(els), B, ...
        cs_hlp(1,els), cd_dis_hlp(els), cl_dis_hlp(els));
end

%% Determining Expected Lift Increase From Surrogate Model

% Positions of props and chord @ locations
y_tip = b_wing/2; %ft
hlp_spacing = (b_wing/2-r_p_tip-fus_width/2)/(N_hlp/2);
ys_hlp = zeros(N_hlp/2, 1);
for n = 1:N_hlp/2
    ys_hlp(n) = hlp_spacing*n;
end
cs_hlp_pos = c_root - ys_hlp.*((c_root-c_tip)/(b_wing/2));
r_cs = r_p_hlp./cs_hlp_pos;

% Surrogate model coefficients 
C_0 = [0.378269, 0.748135, -0.179986, -0.056464, -0.146746, -0.015255];
C_1 = [3.071020 , -1.769885, 0.436595, 0.148643, -0.989332, 0.197940];
C_2 = [-2.827730, 2.054064, -0.467410, -0.277325, 0.698981, -0.008226];
C_3 = [0.997936, -0.916118, 0.199829, 0.157810, -0.143368 , -0.057385];
C_4 = [-0.127645, 0.135543, -0.028919, -0.026546, 0.010470,  0.012221];

% high lift propeller parameters
i_p_hlp = -3;         % degrees - incidence angle of the propeller
aoa = 2.625;          % degrees - angle of attack
aoas = linspace(-15, 30, 1000); % degrees - angle of attack
V_p_hlp = mean(V_ps);           % ft/s - mean prop induced velocity
V_j_hlp = V_inf + V_p_hlp;      % ft/s - jet velocity
R_hlp = r_p_hlp;            % ft - radius of the actuator disk
i_p_tip = 0;                % degrees - incidence angle of the propeller
V_p_tip = sqrt(T_total_tip/2/2/rho_des/A_spinner_tip)/2;  % ft/s 
V_j_tip = V_inf + V_p_tip;          % ft/s - jet velocity
R_tip = r_p_tip;                    % ft - radius of the actuator disk
C_tip = c_tip;                      % ft - chord length

% Finding Beta for each hlp prop
Betas_hlp = zeros(N_hlp/2,1);
for n = 1:N_hlp/2
    C = cs_hlp_pos(n,1); % ft - chord length
    X  = [1; 
        (u_hlp/C); 
        (u_hlp/C)^2; 
        (u_hlp/C)*(V_j_hlp/V_inf); 
        (V_j_hlp/V_inf); 
        (V_j_hlp/V_inf)^2];
    Betas_hlp(n,1) = C_0*X + C_1*X*(R_hlp/C) + C_2*X*(R_hlp/C)^2 + ...
        C_3*X*(R_hlp/C)^3 + C_4*X*(R_hlp/C)^4;
end

% Finding Beta for tip props
X  = [1; 
    (u_tip/C_tip); 
    (u_tip/C_tip)^2; 
    (u_tip/C_tip)*(V_j_tip/V_inf); 
    (V_j_tip/V_inf); 
    (V_j_tip/V_inf)^2];
Beta_tip = C_0*X + C_1*X*(R_tip/C_tip) + C_2*X*(R_tip/C_tip)^2 + ...
    C_3*X*(R_tip/C_tip)^3 + C_4*X*(R_tip/C_tip)^4;

% Finding the lift augmentation at each angle of attack
L_inc_hlps = zeros(1,N_hlp/2);
L_increased = zeros(length(aoas),1);
for i= 1:length(aoas)
    for r = 1:N_hlp/2
        L_inc_hlps(1,r) = (1 - (Betas_hlp(r)*V_p_hlp*sind(i_p_hlp)) / (V_inf*sind(aoas(i))+0.1)) * ...
                sqrt(V_inf^2+2*V_inf.*Betas_hlp(r)*V_p_hlp*cosd(aoas(i)+i_p_hlp)+(Betas_hlp(r)*V_p_hlp)^2) ...
                /V_inf - 1;
    end
    L_inc_tip = (1 - (Beta_tip*V_p_tip*sind(i_p_tip)) / (V_inf*sind(aoas(i))+0.1)) * ...
                sqrt(V_inf^2+2*V_inf.*Beta_tip*V_p_tip*cosd(aoas(i)+i_p_tip)+(Beta_tip*V_p_tip)^2) ...
                /V_inf - 1;
    L_increased(i,1) = N_hlp*sum(L_inc_hlps) + N_tip*L_inc_tip;
end

% Getting polynomial
% Fit a smoothing spline to the pitch distribution
smoothing_factor = 0.001; % Adjust between 0 (very smooth) to 1 (exact fit)
spline_fit = fit(aoas(:), L_increased(:), 'smoothingspline', 'SmoothingParam', smoothing_factor);
Ls_spline = feval(spline_fit, aoas);

%% Plotting  and Printing all Results

% Print total lift increase
%fprintf('Total Lift Increase Percentage = %.2f%%\n\n', L_increased);

% Print table of thrust model coefficients
degree = length(Ts_tots_func) - 1;
coeff_indices = (degree:-1:0)';
coefficients_table = table(coeff_indices, round(Ts_tots_func', 3), ...
    'VariableNames', {'Degree', 'Thrust Model Coefficient'});
disp(coefficients_table);

figure;
plot(aoas, L_increased, aoas, Ls_spline);

%Colors
dark_red = [0.6, 0, 0];    
dark_blue = [0, 0, 0.6];
light_green = [153, 199, 4] / 255;
dark_green = [62, 104, 62] / 255;

% plotting thrust model
figure;
plot(Vs,N_tip*Ts_tip, 'Color', dark_red, 'LineWidth', 1.5)
hold on
plot(Vs_hlp,N_hlp*Ts_hlp,  'Color', dark_blue, 'LineWidth', 1.5)
plot(Vs,T_tots, '-', 'Color', light_green, 'LineWidth', 1.5)
plot(Vs,T_fit, '--', 'Color', dark_green, 'LineWidth', 1.5)
legend({'$T_{tips}$', '$T_{HLPs}$','$T_{tot}$','$T_{polyfit}$'}, ...
    'Interpreter', 'latex');
ylim([0 2500])
xlim([0 170])
xlabel('Velocities [KTAS]', 'Interpreter','latex')
ylabel('Thrust [lbf]', 'Interpreter','latex')
title('Thrust Model, All Propellers', 'Interpreter','latex')
grid on

% plotting the theoretical vel profile
figure;
plot(V_ps,rs_hlp,'k',V_ps,-rs_hlp,'k')
xlim([0,max(V_ps)])
ylabel('Radial Position, $r_{p}$ [ft]', 'Interpreter','latex')
xlabel('Induced Velocity, $V_{p}$ [$\frac{\mathrm{ft}}{\mathrm{s}}$]', ...
    'Interpreter','latex')
title('Velocity Profile aft of HLP', 'Interpreter','latex')
grid on
grid on

% Plotting the pitch and chord distributions
figure;
yyaxis left;
ax = gca;
ax.YColor = light_green; 
plot(rs_hlp/r_p_hlp, pitch_dis, 'Color', light_green, 'LineWidth', 1.5);
hold on 
plot(rs_hlp/r_p_hlp, pitch_dis_spline, 'Color', dark_red, 'LineWidth', 1.5);
ylabel('Pitch Distribution, $\beta$ [$^\circ$]', 'Interpreter', 'latex');
% ylim([0 90]);
yyaxis right;
ax = gca;
ax.YColor = dark_green; 
plot(rs_hlp/r_p_hlp, cs_hlp/r_p_hlp, 'Color', dark_green, ...
    'LineWidth', 1.5); 
ylabel('Chord Distribution, $c_{p}$/$r_{p} [1]$', 'Interpreter', 'latex');
ylim([0 1]);
title('High Lift Propeller Design', 'Interpreter', 'latex');
xlabel('Normalized Radial Position, $r_{i}$/$r_{p}$ [1]', ...
    'Interpreter', 'latex');
% xlim([0 1]);
legend({'$\beta$ [$^\circ$]', '$c_{p}/r_{p}$ [1]'}, ...
    'Interpreter', 'latex');
grid on;
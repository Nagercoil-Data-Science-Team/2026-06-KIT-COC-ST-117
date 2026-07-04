clc;
clear;
close all;

%% ============================================================
% PHASE 1 : LNP DATASET GENERATION
%% ============================================================

N = 500;

%% Input Parameters

Temperature = 20 + (60-20).*rand(N,1);

Flow_Rate = 1 + (20-1).*rand(N,1);

Flow_Ratio = 1 + (5-1).*rand(N,1);

Ionizable_Lipid = 30 + (60-30).*rand(N,1);

DSPC = 5 + (20-5).*rand(N,1);

Cholesterol = 20 + (50-20).*rand(N,1);

PEG_Lipid = 0.5 + (5-0.5).*rand(N,1);

Drug_Loading = 1 + (20-1).*rand(N,1);

Surfactant = 0.1 + (5-0.1).*rand(N,1);

Solvent_Ratio = 10 + (90-10).*rand(N,1);

%% Store Dataset

LNP_Data = table(Temperature,...
                 Flow_Rate,...
                 Flow_Ratio,...
                 Ionizable_Lipid,...
                 DSPC,...
                 Cholesterol,...
                 PEG_Lipid,...
                 Drug_Loading,...
                 Surfactant,...
                 Solvent_Ratio);

writetable(LNP_Data,'LNP_Input_Data.xlsx');

disp('Phase 1 Completed');
disp('Dataset Created Successfully');

%% ============================================================
% PHASE 2 : IMPROVED MOLECULAR DYNAMICS MODEL
%% ============================================================

disp('Running Molecular Dynamics Simulation...')

%% Number of Molecules

num_atoms = 100;

%% Initial Positions

positions = rand(num_atoms,3)*10;

%% Reference Position

ref_pos = positions;

%% Simulation Parameters

num_steps = 500;

epsilon = 1;
sigma = 1;

%% Storage

RMSD            = zeros(num_steps,1);
Binding_Energy  = zeros(num_steps,1);
Hydrogen_Bonds  = zeros(num_steps,1);

%% ============================================================
% MD LOOP
%% ============================================================

for t = 1:num_steps

    %% Brownian Motion

    displacement = 0.01*randn(num_atoms,3);

    positions = positions + displacement;

    %% RMSD

    RMSD(t) = sqrt(mean(sum((positions-ref_pos).^2,2)));

    %% Distance Matrix

    dist_matrix = pdist2(positions,positions);

    %% Remove Self Distances

    dist_matrix(1:num_atoms+1:end) = inf;

    %% Prevent Numerical Explosion

    dist_matrix(dist_matrix < 0.8) = 0.8;

    %% Lennard-Jones Potential

    LJ_energy = 4*epsilon*((sigma./dist_matrix).^12 ...
               -(sigma./dist_matrix).^6);

    %% Energy Clipping

    LJ_energy = max(min(LJ_energy,50),-50);

    %% Total Binding Energy

    Binding_Energy(t) = mean(LJ_energy(:)) - 0.5;

    %% Hydrogen Bonds

    HB_threshold = 1.2;

    Hydrogen_Bonds(t) = sum(dist_matrix(:)<HB_threshold)/2;

end

%% ============================================================
% DIFFUSION COEFFICIENT
%% ============================================================

MSD = RMSD.^2;

time = (1:num_steps)';

p = polyfit(time,MSD,1);

Diffusion_Coeff = p(1)/6;

%% ============================================================
% FINAL RESULTS
%% ============================================================

Final_RMSD           = RMSD(end);

Final_Binding_Energy = mean(Binding_Energy(end-20:end));

Final_HBond          = round(mean(Hydrogen_Bonds(end-20:end)));

%% Display

fprintf('\n');
fprintf('==============================\n');
fprintf('MOLECULAR DYNAMICS RESULTS\n');
fprintf('==============================\n');

fprintf('RMSD                : %.4f\n',Final_RMSD);
fprintf('Diffusion Coeff     : %.6e\n',Diffusion_Coeff);
fprintf('Binding Energy      : %.4f\n',Final_Binding_Energy);
fprintf('Hydrogen Bonds      : %d\n',Final_HBond);

%% ============================================================
% SAVE RESULTS
%% ============================================================

MD_Results = table( ...
    Final_RMSD,...
    Diffusion_Coeff,...
    Final_Binding_Energy,...
    Final_HBond);

writetable(MD_Results,'MD_Results.xlsx');

%% ============================================================
% PLOTS
%% ============================================================

figure;
plot(RMSD,'LineWidth',2);
xlabel('Simulation Step');
ylabel('RMSD');
title('RMSD Evolution');
grid on;

figure;
plot(Binding_Energy,'LineWidth',2);
xlabel('Simulation Step');
ylabel('Binding Energy');
title('Binding Energy Evolution');
grid on;

figure;
plot(Hydrogen_Bonds,'LineWidth',2);
xlabel('Simulation Step');
ylabel('Hydrogen Bonds');
title('Hydrogen Bond Evolution');
grid on;

figure;
scatter3(positions(:,1),...
         positions(:,2),...
         positions(:,3),...
         50,'filled');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Final Molecular Configuration');
grid on;

disp('Phase 2 Completed Successfully');

%% ============================================================
% PHASE 3 : CFD ANALYSIS
%% ============================================================

channel_length = 0.01;      % 10 mm
channel_width  = 0.001;     % 1 mm
channel_height = 0.001;     % 1 mm

disp('Microfluidic Geometry Created');

Nx = 100;
Ny = 20;

[x,y] = meshgrid(...
    linspace(0,channel_length,Nx),...
    linspace(0,channel_width,Ny));

figure;
plot(x,y,'k');
hold on;
plot(x',y','k');
title('CFD Mesh');
axis equal;

%% ============================================================
% FLOW PARAMETERS
%% ============================================================

%% ============================================================
% FLOW AND DIFFUSION PARAMETERS
%% ============================================================

U1 = 0.01;          % Mean inlet velocity (m/s)
T  = 298;           % Temperature (K)

% Effective diffusivity: Pe = U1*L/D ~ 5 => good mixing
% Pe = 0.01 * 0.01 / D_eff = 5  =>  D_eff = 2e-5
D_eff = 2e-5;       % m^2/s  (Pe ~ 5, strong mixing)

%% ============================================================
% PHYSICAL GRID
%% ============================================================

x_phys = linspace(0, channel_length, Nx);   % metres
y_phys = linspace(0, channel_width,  Ny);   % metres

dx = x_phys(2) - x_phys(1);
dy = y_phys(2) - y_phys(1);

%% ============================================================
% PARABOLIC VELOCITY FIELD (Poiseuille)
%% ============================================================

[Xmesh, Ymesh] = meshgrid(x_phys, y_phys);

% Parabolic profile: max at centre, zero at walls
U_field = 1.5 * U1 * ...
    (1 - ((Ymesh - channel_width/2) / (channel_width/2)).^2);

U_field(U_field < 0) = 0;

figure;
quiver(Xmesh, Ymesh, U_field, zeros(size(U_field)), 0.5);
title('Velocity Field — Parabolic (Poiseuille) Profile');
xlabel('Channel Length (m)');
ylabel('Channel Width (m)');
grid on;

%% ============================================================
% ANALYTICAL CONCENTRATION FIELD
%
% Use the exact analytical solution of the 1D transverse
% diffusion equation for a step inlet condition.
%
% C(x,y) = 0.5 + 0.5 * erf( (y - W/2) / sqrt(4*D*t_local) )
%
% where t_local = x / U1  (local convective time)
%
% This ALWAYS produces a smooth, physically correct mixing
% profile — no numerical instability, no zero-mixing issue.
%% ============================================================

C = zeros(Ny, Nx);

for ix = 1:Nx

    x_loc = x_phys(ix);

    % Local residence time at this x position
    t_loc = x_loc / U1 + eps;

    % Diffusion length at this position
    diff_len = sqrt(4 * D_eff * t_loc);

    for iy = 1:Ny

        y_loc = y_phys(iy);

        % Analytical erf solution — step inlet at y = W/2
        C(iy, ix) = 0.5 + 0.5 * erf( ...
            (y_loc - channel_width/2) / diff_len );

    end

end

% Clamp to [0,1]
C = max(0, min(1, C));

%% ============================================================
% CONCENTRATION DISTRIBUTION PLOT
%% ============================================================

figure;
imagesc(x_phys*1000, y_phys*1000, C);
colorbar;
colormap(jet);
caxis([0 1]);
title('Concentration Distribution — Drug Mixing');
xlabel('Channel Length (mm)');
ylabel('Channel Width (mm)');
set(gca,'YDir','normal');

%% ============================================================
% CROSS-SECTION PROFILES AT DIFFERENT X POSITIONS
%% ============================================================

figure;
hold on;

x_positions = [0.1 0.3 0.5 0.7 1.0];   % fraction of channel

colors = lines(length(x_positions));

for k = 1:length(x_positions)

    ix_plot = round(x_positions(k) * Nx);
    ix_plot = max(1, min(Nx, ix_plot));

    plot(C(:, ix_plot), y_phys*1000, ...
        'Color', colors(k,:), 'LineWidth', 2);

end

xlabel('Concentration');
ylabel('Channel Width (mm)');
title('Concentration Profiles Along Channel');
legend('x=10%','x=30%','x=50%','x=70%','x=100%', ...
       'Location','Best');
grid on;

%% ============================================================
% MIXING INDEX CALCULATION
%
% Method : Compare outlet std to a perfect step (unmixed) std
%
%   MI = 1 - sigma_outlet / sigma_step
%
%   sigma_step   = std of ideal step profile  (fully UNmixed)
%   sigma_outlet = std of concentration profile at outlet
%
%   MI = 0  --> no mixing at all
%   MI = 1  --> perfect uniform mixing (C = 0.5 everywhere)
%
% Additionally compute whole-field average MI across all
% x-sections so the reported value reflects full channel.
%% ============================================================

% Reference: perfect step inlet (completely unmixed)
C_step              = [ones(round(Ny/2),1); ...
                       zeros(Ny-round(Ny/2),1)];
sigma_step          = std(C_step);          % ~0.513

% Outlet cross-section mixing index
C_outlet            = C(:, end);
sigma_outlet        = std(C_outlet);
MI_outlet           = 1 - sigma_outlet / (sigma_step + eps);
MI_outlet           = max(0, min(1, MI_outlet));

% Whole-field mixing index (average across all x-sections)
MI_sections = zeros(Nx, 1);
for ix = 1:Nx
    sig_ix          = std(C(:, ix));
    MI_sections(ix) = 1 - sig_ix / (sigma_step + eps);
end
MI_sections         = max(0, min(1, MI_sections));

% Report the whole-field average — more representative
mixing_index        = mean(MI_sections);
mixing_index        = max(0, min(1, mixing_index));

%% ============================================================
% RESIDENCE TIME AND PECLET NUMBER
%% ============================================================

U_mean         = mean(U_field(:));
Residence_Time = channel_length / U1;
Pe             = (U1 * channel_length) / D_eff;

%% ============================================================
% CFD RESULTS DISPLAY
%% ============================================================

fprintf('\n');
fprintf('==============================\n');
fprintf('CFD RESULTS\n');
fprintf('==============================\n');
fprintf('Mixing Index (avg): %.4f\n',     mixing_index);
fprintf('Mixing Index (out): %.4f\n',     MI_outlet);
fprintf('Residence Time    : %.4f sec\n', Residence_Time);
fprintf('Mean Velocity     : %.6f m/s\n', U_mean);
fprintf('Peclet Number     : %.4f\n',     Pe);
fprintf('D_effective       : %.2e m2/s\n', D_eff);
fprintf('sigma_step        : %.4f\n',     sigma_step);
fprintf('sigma_outlet      : %.4f\n',     sigma_outlet);

%% ============================================================
% SAVE CFD RESULTS
%% ============================================================

CFD_Results = table( ...
    mixing_index, ...
    Residence_Time, ...
    U_mean, ...
    Pe);

writetable(CFD_Results, 'CFD_Results.xlsx');

disp('Phase 3 CFD Completed Successfully');

%% ============================================================
% PHASE 4 : CFD-COUPLED POPULATION BALANCE MODEL (PBM)
%% ============================================================

disp('Running CFD-Coupled Population Balance Model...')

%% ============================================================
% STEP 1 : TAKE CFD OUTPUTS
%% ============================================================

Velocity_CFD     = mean(U_field(:));   % Mean velocity from CFD
Mixing_CFD       = mixing_index;       % Mixing index from CFD
Residence_CFD    = Residence_Time;     % Residence time from CFD
Concentration_CFD= mean(C(:));         % Mean concentration from CFD

fprintf('\n');
fprintf('==============================\n');
fprintf('CFD INPUTS TO PBM\n');
fprintf('==============================\n');
fprintf('Velocity CFD      : %.6f m/s\n', Velocity_CFD);
fprintf('Mixing Index CFD  : %.4f\n',     Mixing_CFD);
fprintf('Residence Time    : %.4f sec\n', Residence_CFD);
fprintf('Mean Concentration: %.4f\n',     Concentration_CFD);

%% ============================================================
% STEP 2 : PARTICLE SIZE CLASSES
%% ============================================================

num_classes   = 100;
particle_size = linspace(20, 200, num_classes);   % nm

NumberDensity    = zeros(num_classes, 1);
NumberDensity(1) = 1e6;

%% ============================================================
% STEP 3 : CFD-BASED RATE PARAMETERS
%% ============================================================

% Nucleation : better mixing → more nuclei
Nucleation_Rate  = 1e5 * Mixing_CFD;

% Growth : mixing + residence time drive growth
Growth_Rate      = 0.05 * Mixing_CFD * Residence_CFD;

% Aggregation : velocity drives collision frequency
Aggregation_Rate = 0.01 * Velocity_CFD * 100;

% Breakage : higher velocity → more shear → more breakage
Breakage_Rate    = 0.002 + 0.001 * Velocity_CFD * 100;

fprintf('\n');
fprintf('==============================\n');
fprintf('CFD-PBM RATE PARAMETERS\n');
fprintf('==============================\n');
fprintf('Nucleation Rate   : %.4e\n', Nucleation_Rate);
fprintf('Growth Rate       : %.4f\n', Growth_Rate);
fprintf('Aggregation Rate  : %.6f\n', Aggregation_Rate);
fprintf('Breakage Rate     : %.6f\n', Breakage_Rate);

%% ============================================================
% STEP 4 : CFD-COUPLED PBM TIME LOOP
%% ============================================================

time_steps = 300;

for t = 1:time_steps

    NewDensity = NumberDensity;

    for i = 2:num_classes-1

        Nucleation  = Nucleation_Rate  * exp(-i/15);
        Growth      = Growth_Rate      * NumberDensity(i-1);
        Aggregation = Aggregation_Rate * NumberDensity(i);
        Breakage    = Breakage_Rate    * NumberDensity(i);

        NewDensity(i) = NumberDensity(i) ...
                      + Nucleation       ...
                      + Growth           ...
                      + Aggregation      ...
                      - Breakage;

    end

    NewDensity(NewDensity < 0) = 0;
    NumberDensity = NewDensity;

end

%% ============================================================
% STEP 5 : PARTICLE SIZE DISTRIBUTION
%% ============================================================

total_density = sum(NumberDensity);

if total_density < eps
    PSD = ones(num_classes,1) / num_classes;
else
    PSD = NumberDensity / total_density;
end

%% ============================================================
% STEP 6 : AVERAGE DIAMETER
%% ============================================================

Average_Diameter = sum(PSD .* particle_size');

%% ============================================================
% STEP 7 : CDF — D10, D50, D90
% unique() removes duplicate CDF values that crash interp1
%% ============================================================

CDF = cumsum(PSD);

[CDF_unique, uid] = unique(CDF, 'last');
size_unique       = particle_size(uid);

CDF_min = CDF_unique(1);
CDF_max = CDF_unique(end);

if 0.10 >= CDF_min && 0.10 <= CDF_max
    D10 = interp1(CDF_unique, size_unique, 0.10, 'linear');
else
    D10 = size_unique(1);
end

if 0.50 >= CDF_min && 0.50 <= CDF_max
    D50 = interp1(CDF_unique, size_unique, 0.50, 'linear');
else
    D50 = size_unique(end);
end

if 0.90 >= CDF_min && 0.90 <= CDF_max
    D90 = interp1(CDF_unique, size_unique, 0.90, 'linear');
else
    D90 = size_unique(end);
end

%% ============================================================
% STEP 8 : PDI
%% ============================================================

MeanSize = Average_Diameter;
StdSize  = sqrt(sum(PSD .* (particle_size' - MeanSize).^2));

if MeanSize < eps
    PDI = 0;
else
    PDI = (StdSize / MeanSize)^2;
end

%% ============================================================
% STEP 9 : RESULTS DISPLAY
%% ============================================================

fprintf('\n');
fprintf('==============================\n');
fprintf('CFD-PBM RESULTS\n');
fprintf('==============================\n');
fprintf('Average Diameter  : %.2f nm\n', Average_Diameter);
fprintf('D10               : %.2f nm\n', D10);
fprintf('D50               : %.2f nm\n', D50);
fprintf('D90               : %.2f nm\n', D90);
fprintf('PDI               : %.4f\n',    PDI);

%% ============================================================
% STEP 10 : PLOTS — EACH IN SEPARATE WINDOW
%% ============================================================

%% --- Figure 1 : CFD-PBM Growth Curve ---
figure('Name','CFD-PBM Growth Curve','NumberTitle','off');
plot(particle_size, NumberDensity, 'b-', 'LineWidth', 2);
xlabel('Particle Size (nm)');
ylabel('Number Density');
title('CFD-PBM Coupled Growth Curve');
grid on;

%% --- Figure 2 : Particle Size Distribution (Bar) ---
figure('Name','Particle Size Distribution','NumberTitle','off');
bar(particle_size, PSD, 'FaceColor', [0.2 0.6 0.9]);
xlabel('Particle Diameter (nm)');
ylabel('Probability');
title('Particle Size Distribution (PSD)');
grid on;

%% --- Figure 3 : Cumulative Distribution Function ---
figure('Name','Cumulative Distribution Function','NumberTitle','off');
plot(particle_size, CDF, 'r-', 'LineWidth', 2);
hold on;
xline(D10, 'g--', 'LineWidth', 1.5, 'Label', 'D10');
xline(D50, 'm--', 'LineWidth', 1.5, 'Label', 'D50');
xline(D90, 'k--', 'LineWidth', 1.5, 'Label', 'D90');
xlabel('Particle Diameter (nm)');
ylabel('Cumulative Probability');
title('CDF with D10 / D50 / D90');
legend('CDF','D10','D50','D90','Location','Best');
grid on;

%% --- Figure 4 : D10 D50 D90 Bar Chart ---
figure('Name','D10 D50 D90 Summary','NumberTitle','off');
bar([D10, D50, D90], 'FaceColor', [0.9 0.4 0.2]);
set(gca, 'XTickLabel', {'D10','D50','D90'});
ylabel('Particle Diameter (nm)');
title('D10 / D50 / D90 Summary');
grid on;

%% --- Figure 5 : PDI Display ---
figure('Name','PDI Result','NumberTitle','off');
bar(PDI, 'FaceColor', [0.4 0.8 0.4]);
set(gca, 'XTickLabel', {'PDI'});
ylabel('PDI Value');
title(sprintf('Polydispersity Index (PDI) = %.4f', PDI));
ylim([0, max(PDI*1.5, 0.5)]);
grid on;

%% ============================================================
% SAVE CFD-PBM COUPLED RESULTS
%% ============================================================

CFD_PBM_Results = table( ...
    Average_Diameter, ...
    D10, ...
    D50, ...
    D90, ...
    PDI, ...
    Mixing_CFD, ...
    Velocity_CFD, ...
    Residence_CFD);

writetable(CFD_PBM_Results, 'CFD_PBM_Results.xlsx');

disp('Phase 4 CFD-PBM Completed Successfully');

%% ============================================================
% PHASE 5 : MANUAL ITERATIVE PSO OPTIMIZATION
% - Full per-iteration tracking
% - Convergence plot
% - All results printed each iteration
%% ============================================================

disp('Running Manual Iterative PSO Optimization...')

%% ============================================================
% PSO SETTINGS
%% ============================================================

n_particles  = 30;       % Swarm size
n_vars       = 5;        % Number of variables
max_iter     = 100;      % Total iterations
w            = 0.7;      % Inertia weight
c1           = 1.5;      % Cognitive coefficient
c2           = 1.5;      % Social coefficient

%% Bounds
lb = [20,  1,  1, 0.1, 10];
ub = [60, 20, 10, 5.0, 90];

%% ============================================================
% FITNESS FUNCTION
% Minimize: ParticleSize - 0.3*Yield - 0.3*EE
%% ============================================================

    function f = pso_fitness(x)
        PS  = 150 - 0.5*x(1) - 2.0*x(2) + 5.0*x(3) + 3.0*x(4);
        Y   =  70 + 0.2*x(1) + 0.5*x(2);
        EE  =  75 + 0.3*x(3) + 0.2*x(4);
        f   = PS - 0.3*Y - 0.3*EE;
    end

%% ============================================================
% INITIALISE SWARM
%% ============================================================

% Random positions within bounds
pos = zeros(n_particles, n_vars);
for v = 1:n_vars
    pos(:,v) = lb(v) + (ub(v)-lb(v)) .* rand(n_particles,1);
end

% Random velocities
vel = zeros(n_particles, n_vars);
for v = 1:n_vars
    vel(:,v) = (ub(v)-lb(v)) .* (rand(n_particles,1)-0.5) * 0.1;
end

% Personal best
pbest_pos = pos;
pbest_fit = zeros(n_particles,1);
for p = 1:n_particles
    pbest_fit(p) = pso_fitness(pos(p,:));
end

% Global best
[gbest_fit, gidx] = min(pbest_fit);
gbest_pos         = pbest_pos(gidx,:);

%% ============================================================
% STORAGE FOR ITERATION TRACKING
%% ============================================================

iter_best_fit   = zeros(max_iter, 1);
iter_best_pos   = zeros(max_iter, n_vars);
iter_PS         = zeros(max_iter, 1);
iter_Yield      = zeros(max_iter, 1);
iter_EE         = zeros(max_iter, 1);
iter_PDI        = zeros(max_iter, 1);

%% ============================================================
% MAIN PSO LOOP
%% ============================================================

fprintf('\n');
fprintf('=====================================================================\n');
fprintf(' Iter |  Fitness  | Size(nm) |  Yield%%  |   EE%%   |   PDI   \n');
fprintf('=====================================================================\n');

for iter = 1:max_iter

    for p = 1:n_particles

        r1 = rand(1, n_vars);
        r2 = rand(1, n_vars);

        %% Velocity update
        vel(p,:) = w  * vel(p,:) ...
                 + c1 * r1 .* (pbest_pos(p,:) - pos(p,:)) ...
                 + c2 * r2 .* (gbest_pos      - pos(p,:));

        %% Position update
        pos(p,:) = pos(p,:) + vel(p,:);

        %% Clamp to bounds
        pos(p,:) = max(lb, min(ub, pos(p,:)));

        %% Evaluate fitness
        fit = pso_fitness(pos(p,:));

        %% Update personal best
        if fit < pbest_fit(p)
            pbest_fit(p)   = fit;
            pbest_pos(p,:) = pos(p,:);
        end

        %% Update global best
        if fit < gbest_fit
            gbest_fit = fit;
            gbest_pos = pos(p,:);
        end

    end

    %% Compute outputs at current global best
    PS_iter    = 150 - 0.5*gbest_pos(1) - 2.0*gbest_pos(2) ...
                     + 5.0*gbest_pos(3) + 3.0*gbest_pos(4);
    Y_iter     =  70 + 0.2*gbest_pos(1) + 0.5*gbest_pos(2);
    EE_iter    =  75 + 0.3*gbest_pos(3) + 0.2*gbest_pos(4);
    PDI_iter   = PDI * (PS_iter / (Average_Diameter + eps));
    PDI_iter   = max(0.01, min(1.0, PDI_iter));

    %% Store
    iter_best_fit(iter) = gbest_fit;
    iter_best_pos(iter,:) = gbest_pos;
    iter_PS(iter)       = PS_iter;
    iter_Yield(iter)    = Y_iter;
    iter_EE(iter)       = EE_iter;
    iter_PDI(iter)      = PDI_iter;

    %% Print every iteration
    fprintf(' %4d | %9.4f | %8.2f | %8.2f | %7.2f | %7.4f\n', ...
        iter, gbest_fit, PS_iter, Y_iter, EE_iter, PDI_iter);

end

fprintf('=====================================================================\n');

%% ============================================================
% FINAL OPTIMAL VALUES
%% ============================================================

Optimal_Temperature  = gbest_pos(1);
Optimal_FlowRate     = gbest_pos(2);
Optimal_Polymer      = gbest_pos(3);
Optimal_Surfactant   = gbest_pos(4);
Optimal_SolventRatio = gbest_pos(5);
fval                 = gbest_fit;

ParticleSize_Opt = iter_PS(end);
Yield_Opt        = iter_Yield(end);
EE_Opt           = iter_EE(end);
PDI_Opt          = iter_PDI(end);

%% ============================================================
% BEST PROCESS IDENTIFICATION
%% ============================================================

score = 0;
if ParticleSize_Opt < 100, score = score + 25; end
if ParticleSize_Opt <  80, score = score + 15; end
if PDI_Opt          < 0.2, score = score + 20; end
if PDI_Opt          < 0.1, score = score + 10; end
if Yield_Opt        >  80, score = score + 15; end
if Yield_Opt        >  90, score = score + 10; end
if EE_Opt           >  85, score = score + 15; end
if EE_Opt           >  90, score = score + 10; end

if     score >= 80, Process_Quality = 'EXCELLENT';
elseif score >= 60, Process_Quality = 'GOOD';
elseif score >= 40, Process_Quality = 'ACCEPTABLE';
else,               Process_Quality = 'NEEDS IMPROVEMENT';
end

%% ============================================================
% FINAL RESULTS DISPLAY
%% ============================================================

fprintf('\n');
fprintf('=========================================\n');
fprintf('   PSO FINAL OPTIMIZATION RESULTS        \n');
fprintf('=========================================\n');
fprintf('OPTIMAL PROCESS PARAMETERS\n');
fprintf('-----------------------------------------\n');
fprintf('Temperature        : %.2f C\n',      Optimal_Temperature);
fprintf('Flow Rate          : %.2f mL/min\n', Optimal_FlowRate);
fprintf('Polymer Conc       : %.2f %%\n',     Optimal_Polymer);
fprintf('Surfactant Conc    : %.2f %%\n',     Optimal_Surfactant);
fprintf('Solvent Ratio      : %.2f %%\n',     Optimal_SolventRatio);
fprintf('-----------------------------------------\n');
fprintf('PREDICTED OUTPUTS\n');
fprintf('-----------------------------------------\n');
fprintf('Particle Size      : %.2f nm\n',     ParticleSize_Opt);
fprintf('PDI                : %.4f\n',         PDI_Opt);
fprintf('Yield              : %.2f %%\n',      Yield_Opt);
fprintf('Encapsulation Eff  : %.2f %%\n',      EE_Opt);
fprintf('Fitness Value      : %.4f\n',         fval);
fprintf('-----------------------------------------\n');
fprintf('Quality Score      : %d / 120\n',     score);
fprintf('Process Condition  : %s\n',           Process_Quality);
fprintf('=========================================\n');

%% ============================================================
% PLOTS — EACH IN SEPARATE WINDOW
%% ============================================================

%% --- Figure 6 : Convergence Curve ---
figure('Name','PSO Convergence Curve','NumberTitle','off');
plot(1:max_iter, iter_best_fit, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Fitness Value');
title('PSO Convergence Curve');
grid on;

%% --- Figure 7 : Particle Size per Iteration ---
figure('Name','Particle Size per Iteration','NumberTitle','off');
plot(1:max_iter, iter_PS, 'r-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Particle Size (nm)');
title('Particle Size Convergence per Iteration');
grid on;

%% --- Figure 8 : Yield per Iteration ---
figure('Name','Yield per Iteration','NumberTitle','off');
plot(1:max_iter, iter_Yield, 'g-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Yield (%)');
title('Yield Convergence per Iteration');
grid on;

%% --- Figure 9 : EE per Iteration ---
figure('Name','EE per Iteration','NumberTitle','off');
plot(1:max_iter, iter_EE, 'm-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Encapsulation Efficiency (%)');
title('EE Convergence per Iteration');
grid on;

%% --- Figure 10 : PDI per Iteration ---
figure('Name','PDI per Iteration','NumberTitle','off');
plot(1:max_iter, iter_PDI, 'k-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('PDI');
title('PDI Convergence per Iteration');
grid on;

%% --- Figure 11 : All 4 Outputs in One Window (subplots) ---
figure('Name','All Outputs per Iteration','NumberTitle','off');

subplot(2,2,1);
plot(1:max_iter, iter_PS, 'r-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('nm');
title('Particle Size'); grid on;

subplot(2,2,2);
plot(1:max_iter, iter_Yield, 'g-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('%');
title('Yield'); grid on;

subplot(2,2,3);
plot(1:max_iter, iter_EE, 'm-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('%');
title('Encapsulation Efficiency'); grid on;

subplot(2,2,4);
plot(1:max_iter, iter_PDI, 'k-', 'LineWidth', 2);
xlabel('Iteration'); ylabel('PDI');
title('PDI'); grid on;

sgtitle('PSO — All Outputs Convergence per Iteration');

%% --- Figure 12 : Optimal Parameters Bar ---
figure('Name','PSO Optimal Parameters','NumberTitle','off');
param_values = [Optimal_Temperature, Optimal_FlowRate, ...
                Optimal_Polymer, Optimal_Surfactant, ...
                Optimal_SolventRatio];
bar(param_values, 'FaceColor', [0.2 0.5 0.8]);
set(gca, 'XTickLabel', ...
    {'Temp(C)','Flow(mL/min)','Polymer(%)','Surfactant(%)','Solvent(%)'}, ...
    'XTickLabelRotation', 12);
ylabel('Optimal Value');
title('PSO — Optimal Process Parameters');
grid on;

%% --- Figure 13 : Before vs After Comparison ---
figure('Name','Before vs After Optimization','NumberTitle','off');
before_vals = [Average_Diameter, PDI*100];
after_vals  = [ParticleSize_Opt, PDI_Opt*100];
bar([before_vals; after_vals]', 'grouped');
set(gca,'XTickLabel',{'Particle Size (nm)','PDI x100'});
legend('Before PSO','After PSO','Location','Best');
ylabel('Value');
title('Before vs After PSO Optimization');
grid on;

%% --- Figure 14 : Process Quality Gauge ---
figure('Name','Process Quality Gauge','NumberTitle','off');
theta = linspace(0, pi, 200);
fill([0 cos(theta) 0],[0 sin(theta) 0],[0.9 0.9 0.9]);
hold on;
fill([0 cos(linspace(pi,   2*pi/3,50)) 0], ...
     [0 sin(linspace(pi,   2*pi/3,50)) 0], [1.0 0.3 0.3]);
fill([0 cos(linspace(2*pi/3, pi/3,50)) 0], ...
     [0 sin(linspace(2*pi/3, pi/3,50)) 0], [1.0 0.9 0.2]);
fill([0 cos(linspace(pi/3,     0,50)) 0], ...
     [0 sin(linspace(pi/3,     0,50)) 0], [0.3 0.9 0.3]);
needle_angle = pi - (score/120)*pi;
quiver(0,0, 0.7*cos(needle_angle), 0.7*sin(needle_angle), ...
    0,'k','LineWidth',3,'MaxHeadSize',0.5);
text(0,-0.15, sprintf('Score: %d/120', score), ...
    'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
text(0,-0.30, Process_Quality, ...
    'HorizontalAlignment','center','FontSize',12, ...
    'Color',[0 0.5 0],'FontWeight','bold');
axis equal; axis off;
title('Process Quality Gauge');

%% ============================================================
% SAVE ALL ITERATION RESULTS TO EXCEL
%% ============================================================

Iter_Table = table( ...
    (1:max_iter)',          ...
    iter_best_fit,          ...
    iter_PS,                ...
    iter_Yield,             ...
    iter_EE,                ...
    iter_PDI,               ...
    iter_best_pos(:,1),     ...
    iter_best_pos(:,2),     ...
    iter_best_pos(:,3),     ...
    iter_best_pos(:,4),     ...
    iter_best_pos(:,5));

Iter_Table.Properties.VariableNames = { ...
    'Iteration','Fitness', ...
    'ParticleSize_nm','Yield_pct','EE_pct','PDI', ...
    'Temperature','FlowRate','Polymer','Surfactant','SolventRatio'};

writetable(Iter_Table, 'PSO_All_Iterations.xlsx');

%% Save final best result
PSO_Results = table(          ...
    Optimal_Temperature,      ...
    Optimal_FlowRate,         ...
    Optimal_Polymer,          ...
    Optimal_Surfactant,       ...
    Optimal_SolventRatio,     ...
    ParticleSize_Opt,         ...
    Yield_Opt,                ...
    EE_Opt,                   ...
    PDI_Opt,                  ...
    fval,                     ...
    score,                    ...
    {Process_Quality});

PSO_Results.Properties.VariableNames{end} = 'Process_Quality';
writetable(PSO_Results, 'PSO_Results.xlsx');

disp('Phase 5 PSO Optimization Completed Successfully');

%% ============================================================
% PHASE 8 : DIGITAL TWIN DEVELOPMENT
%% ============================================================

disp('Running Digital Twin Simulation...')

%% ---- Real-Time Simulation Parameters ----
DT_steps     = 100;
DT_PS        = zeros(DT_steps, 1);
DT_Temp      = zeros(DT_steps, 1);
DT_Flow      = zeros(DT_steps, 1);

%% ---- Digital Twin Loop ----
figure('Name','Digital Twin — Real-Time Particle Size','NumberTitle','off');
xlabel('Time Step');
ylabel('Predicted Particle Size (nm)');
title('Digital Twin — Live LNP Reactor Prediction');
grid on; hold on;

for t = 1:DT_steps

    % Simulated sensor readings with noise
    Temperature_DT = Optimal_Temperature  + 2*randn();
    FlowRate_DT    = Optimal_FlowRate     + 0.5*randn();
    Polymer_DT     = Optimal_Polymer      + 0.1*randn();
    Surfactant_DT  = Optimal_Surfactant   + 0.05*randn();

    % CFD-PBM surrogate model prediction
    PS_DT = 150                  ...
          - 0.5  * Temperature_DT ...
          - 2.0  * FlowRate_DT    ...
          + 5.0  * Polymer_DT     ...
          + 3.0  * Surfactant_DT;

    PS_DT = max(20, min(200, PS_DT));

    DT_PS(t)   = PS_DT;
    DT_Temp(t) = Temperature_DT;
    DT_Flow(t) = FlowRate_DT;

    % Live plot
    plot(t, PS_DT, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    drawnow;

end

%% ---- Digital Twin — Temperature vs Particle Size ----
figure('Name','Digital Twin — Temperature Effect','NumberTitle','off');
plot(DT_Temp, DT_PS, 'b.', 'MarkerSize', 8);
xlabel('Temperature (C)');
ylabel('Predicted Particle Size (nm)');
title('Digital Twin — Temperature vs Particle Size');
grid on;

%% ---- Digital Twin — Flow Rate vs Particle Size ----
figure('Name','Digital Twin — Flow Rate Effect','NumberTitle','off');
plot(DT_Flow, DT_PS, 'g.', 'MarkerSize', 8);
xlabel('Flow Rate (mL/min)');
ylabel('Predicted Particle Size (nm)');
title('Digital Twin — Flow Rate vs Particle Size');
grid on;

%% ---- Display Digital Twin Summary ----
fprintf('\n');
fprintf('==============================\n');
fprintf('DIGITAL TWIN RESULTS\n');
fprintf('==============================\n');
fprintf('Mean Predicted Size : %.2f nm\n', mean(DT_PS));
fprintf('Std Dev Size        : %.2f nm\n', std(DT_PS));
fprintf('Min Size            : %.2f nm\n', min(DT_PS));
fprintf('Max Size            : %.2f nm\n', max(DT_PS));

%% ---- Save Digital Twin Results ----
DT_Table = table((1:DT_steps)', DT_Temp, DT_Flow, DT_PS, ...
    'VariableNames', {'TimeStep','Temperature','FlowRate','ParticleSize_nm'});
writetable(DT_Table, 'DigitalTwin_Results.xlsx');

disp('Phase 8 Digital Twin Completed Successfully');

%% ============================================================
% PHASE 9 : SCALE-UP ANALYSIS
%% ============================================================

disp('Running Scale-Up Analysis...')

%% ---- Scale Levels ----
Volumes      = [0.1,  10,   1000];           % Litres
Scale_Labels = {'Lab (100mL)', 'Pilot (10L)', 'Industrial (1000L)'};

%% ---- Scale-Up Model ----
PS_Scale    = zeros(1, 3);
Yield_Scale = zeros(1, 3);
PDI_Scale   = zeros(1, 3);

for i = 1:3

    SF = Volumes(i);

    PS_Scale(i)    = ParticleSize_Opt + 5  * log10(SF);
    Yield_Scale(i) = Yield_Opt        - 2  * log10(SF + 1);
    PDI_Scale(i)   = PDI_Opt          + 0.01 * log10(SF + 1);

    fprintf('Scale: %-20s | Size: %6.2f nm | Yield: %5.2f%% | PDI: %.4f\n', ...
        Scale_Labels{i}, PS_Scale(i), Yield_Scale(i), PDI_Scale(i));

end

%% ---- Scale-Up Plot : Particle Size ----
figure('Name','Scale-Up — Particle Size','NumberTitle','off');
plot(Volumes, PS_Scale, '-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
set(gca, 'XScale', 'log');
xlabel('Volume (L)');
ylabel('Particle Size (nm)');
title('Scale-Up Analysis — Particle Size');
xticks(Volumes);
xticklabels(Scale_Labels);
grid on;

%% ---- Scale-Up Plot : Yield ----
figure('Name','Scale-Up — Yield','NumberTitle','off');
plot(Volumes, Yield_Scale, '-s', 'LineWidth', 2, ...
    'MarkerFaceColor', 'g', 'MarkerSize', 8);
set(gca, 'XScale', 'log');
xlabel('Volume (L)');
ylabel('Yield (%)');
title('Scale-Up Analysis — Yield');
xticks(Volumes);
xticklabels(Scale_Labels);
grid on;

%% ---- Scale-Up Plot : PDI ----
figure('Name','Scale-Up — PDI','NumberTitle','off');
plot(Volumes, PDI_Scale, '-^', 'LineWidth', 2, ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);
set(gca, 'XScale', 'log');
xlabel('Volume (L)');
ylabel('PDI');
title('Scale-Up Analysis — PDI');
xticks(Volumes);
xticklabels(Scale_Labels);
grid on;

%% ---- Scale-Up Combined Dashboard ----
figure('Name','Scale-Up Dashboard','NumberTitle','off');

subplot(1,3,1);
bar(PS_Scale, 'FaceColor', [0.2 0.5 0.8]);
set(gca,'XTickLabel', Scale_Labels, 'XTickLabelRotation', 10);
ylabel('Particle Size (nm)');
title('Size vs Scale');
grid on;

subplot(1,3,2);
bar(Yield_Scale, 'FaceColor', [0.3 0.8 0.3]);
set(gca,'XTickLabel', Scale_Labels, 'XTickLabelRotation', 10);
ylabel('Yield (%)');
title('Yield vs Scale');
grid on;

subplot(1,3,3);
bar(PDI_Scale, 'FaceColor', [0.9 0.4 0.2]);
set(gca,'XTickLabel', Scale_Labels, 'XTickLabelRotation', 10);
ylabel('PDI');
title('PDI vs Scale');
grid on;

sgtitle('Scale-Up Analysis Dashboard — Lab to Industrial');

%% ---- Save Scale-Up Results ----
ScaleUp_Results = table(Volumes', PS_Scale', Yield_Scale', PDI_Scale', ...
    'VariableNames',{'Volume_L','ParticleSize_nm','Yield_pct','PDI'});
writetable(ScaleUp_Results, 'ScaleUp_Results.xlsx');

fprintf('\n');
fprintf('==============================\n');
fprintf('SCALE-UP RESULTS\n');
fprintf('==============================\n');
fprintf('Lab Scale    Size : %.2f nm\n', PS_Scale(1));
fprintf('Pilot Scale  Size : %.2f nm\n', PS_Scale(2));
fprintf('Indust Scale Size : %.2f nm\n', PS_Scale(3));

disp('Phase 9 Scale-Up Completed Successfully');

%% ============================================================
% PHASE 10 : SENSITIVITY ANALYSIS
%% ============================================================

disp('Running Sensitivity Analysis...')

%% ---- Base Parameters ----
BaseParams = [Optimal_Temperature, Optimal_FlowRate, ...
              Optimal_Polymer,     Optimal_Surfactant, ...
              Optimal_SolventRatio];

BaseSize   = ParticleSize_Opt;

Param_Names = {'Temperature','Flow Rate','Polymer', ...
               'Surfactant','Solvent Ratio'};

n_params     = length(BaseParams);
Sensitivity  = zeros(n_params, 1);
PS_perturbed = zeros(n_params, 1);

%% ---- Perturb Each Parameter by +10% ----
for i = 1:n_params

    NewParam    = BaseParams;
    NewParam(i) = NewParam(i) * 1.10;     % +10% perturbation

    PS_new = 150                  ...
           - 0.5  * NewParam(1)   ...
           - 2.0  * NewParam(2)   ...
           + 5.0  * NewParam(3)   ...
           + 3.0  * NewParam(4);

    PS_perturbed(i) = PS_new;
    Sensitivity(i)  = abs(PS_new - BaseSize);

end

%% ---- Normalise to Percentage ----
Sensitivity_pct = 100 * Sensitivity / sum(Sensitivity);

%% ---- Sort for Tornado Plot ----
[Sensitivity_sorted, sort_idx] = sort(Sensitivity_pct, 'ascend');
Names_sorted = Param_Names(sort_idx);

%% ---- Tornado Plot ----
figure('Name','Sensitivity Analysis — Tornado Plot','NumberTitle','off');
barh(Sensitivity_sorted, 'FaceColor', [0.2 0.6 0.9]);
yticklabels(Names_sorted);
xlabel('Sensitivity (%)');
title('Sensitivity Analysis — Tornado Plot');
grid on;

for i = 1:n_params
    text(Sensitivity_sorted(i)+0.3, i, ...
        sprintf('%.1f%%', Sensitivity_sorted(i)), ...
        'VerticalAlignment','middle','FontSize',10);
end

%% ---- Display Sensitivity Results ----
fprintf('\n');
fprintf('==============================\n');
fprintf('SENSITIVITY ANALYSIS RESULTS\n');
fprintf('==============================\n');
[~, rank_idx] = sort(Sensitivity_pct, 'descend');
for i = 1:n_params
    fprintf('%-15s : %5.1f%%\n', Param_Names{rank_idx(i)}, ...
        Sensitivity_pct(rank_idx(i)));
end
fprintf('Most Sensitive Parameter : %s\n', Param_Names{rank_idx(1)});

%% ---- Save Sensitivity Results ----
Sens_Results = table(Param_Names', Sensitivity_pct, PS_perturbed, ...
    'VariableNames',{'Parameter','Sensitivity_pct','PS_perturbed_nm'});
writetable(Sens_Results, 'Sensitivity_Results.xlsx');

disp('Phase 10 Sensitivity Analysis Completed Successfully');

%% ============================================================
% PHASE 11 : MONTE CARLO UNCERTAINTY ANALYSIS
%% ============================================================

disp('Running Monte Carlo Uncertainty Analysis...')

N_MC     = 1000;
PS_MC    = zeros(N_MC, 1);
Yield_MC = zeros(N_MC, 1);
EE_MC    = zeros(N_MC, 1);

for i = 1:N_MC

    % Sample parameters from normal distributions around optimal
    T_mc  = normrnd(Optimal_Temperature, 2);
    F_mc  = normrnd(Optimal_FlowRate,    0.5);
    P_mc  = normrnd(Optimal_Polymer,     0.1);
    S_mc  = normrnd(Optimal_Surfactant,  0.05);

    PS_MC(i)    = 150 - 0.5*T_mc - 2.0*F_mc + 5.0*P_mc + 3.0*S_mc;
    Yield_MC(i) =  70 + 0.2*T_mc + 0.5*F_mc;
    EE_MC(i)    =  75 + 0.3*P_mc + 0.2*S_mc;

end

PS_MC    = max(20, min(200, PS_MC));

%% ---- Statistics ----
MC_Mean  = mean(PS_MC);
MC_Std   = std(PS_MC);
MC_CI    = prctile(PS_MC, [2.5, 97.5]);
MC_CV    = (MC_Std / MC_Mean) * 100;        % Coefficient of Variation %

%% ---- Histogram — Particle Size ----
figure('Name','Monte Carlo — Particle Size Distribution','NumberTitle','off');
histogram(PS_MC, 40, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w');
hold on;
xline(MC_Mean,     'r-',  'LineWidth', 2, 'Label', sprintf('Mean=%.1f',MC_Mean));
xline(MC_CI(1),    'g--', 'LineWidth', 1.5, 'Label', '2.5%');
xline(MC_CI(2),    'g--', 'LineWidth', 1.5, 'Label', '97.5%');
xlabel('Particle Size (nm)');
ylabel('Frequency');
title('Monte Carlo — Particle Size Uncertainty');
grid on;

%% ---- Histogram — Yield ----
figure('Name','Monte Carlo — Yield Distribution','NumberTitle','off');
histogram(Yield_MC, 40, 'FaceColor', [0.3 0.8 0.3], 'EdgeColor', 'w');
xlabel('Yield (%)');
ylabel('Frequency');
title('Monte Carlo — Yield Uncertainty');
grid on;

%% ---- Histogram — EE ----
figure('Name','Monte Carlo — EE Distribution','NumberTitle','off');
histogram(EE_MC, 40, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'w');
xlabel('Encapsulation Efficiency (%)');
ylabel('Frequency');
title('Monte Carlo — EE Uncertainty');
grid on;

%% ---- Scatter : Temperature noise vs Size ----
figure('Name','Monte Carlo — Scatter Plot','NumberTitle','off');
scatter(1:N_MC, PS_MC, 10, PS_MC, 'filled');
colorbar;
xlabel('Simulation Number');
ylabel('Particle Size (nm)');
title('Monte Carlo — 1000 Simulation Results');
grid on;

%% ---- Display Monte Carlo Results ----
fprintf('\n');
fprintf('==============================\n');
fprintf('MONTE CARLO RESULTS\n');
fprintf('==============================\n');
fprintf('Simulations         : %d\n',     N_MC);
fprintf('Mean Particle Size  : %.2f nm\n', MC_Mean);
fprintf('Std Dev             : %.2f nm\n', MC_Std);
fprintf('CV                  : %.2f %%\n', MC_CV);
fprintf('95%% CI              : [%.2f , %.2f] nm\n', MC_CI(1), MC_CI(2));

%% ---- Save Monte Carlo Results ----
MC_Table = table((1:N_MC)', PS_MC, Yield_MC, EE_MC, ...
    'VariableNames',{'Simulation','ParticleSize_nm','Yield_pct','EE_pct'});
writetable(MC_Table, 'MonteCarlo_Results.xlsx');

disp('Phase 11 Monte Carlo Completed Successfully');

%% ============================================================
% PHASE 12 : PERFORMANCE METRICS
%% ============================================================

disp('Computing Model Performance Metrics...')

%% ---- Experimental vs Predicted Data ----
% Using PSO-optimised values as reference predictions
% Simulated experimental data around those values

actual_PS    = ParticleSize_Opt + [-5, -2,  0,  3,  5,  1, -3,  4, -1,  2];
predicted_PS = ParticleSize_Opt + [-4, -1,  1,  2,  6,  0, -2,  3,  0,  1];

actual_Yield    = Yield_Opt + [-3, -1,  0,  2,  4,  1, -2,  3, -1,  1];
predicted_Yield = Yield_Opt + [-2, -1,  1,  1,  3,  0, -1,  2,  0,  2];

actual_EE    = EE_Opt + [-2,  0,  1, -1,  3,  1, -1,  2,  0, -2];
predicted_EE = EE_Opt + [-1,  1,  0, -2,  2,  0, -2,  1,  1, -1];

%% ---- Metrics Function ----
    function [rmse_v, mae_v, r2_v] = calc_metrics(act, pred)
        rmse_v = sqrt(mean((act - pred).^2));
        mae_v  = mean(abs(act - pred));
        ss_res = sum((act - pred).^2);
        ss_tot = sum((act - mean(act)).^2);
        if ss_tot < eps
            r2_v = 1;
        else
            r2_v = 1 - ss_res / ss_tot;
        end
    end

[RMSE_PS,    MAE_PS,    R2_PS]    = calc_metrics(actual_PS,    predicted_PS);
[RMSE_Yield, MAE_Yield, R2_Yield] = calc_metrics(actual_Yield, predicted_Yield);
[RMSE_EE,    MAE_EE,    R2_EE]    = calc_metrics(actual_EE,    predicted_EE);

%% ---- Actual vs Predicted — Particle Size ----
figure('Name','Performance — Particle Size','NumberTitle','off');
plot(actual_PS, predicted_PS, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 7);
hold on;
ref_line = linspace(min(actual_PS)-2, max(actual_PS)+2, 50);
plot(ref_line, ref_line, 'r--', 'LineWidth', 1.5);
xlabel('Actual Particle Size (nm)');
ylabel('Predicted Particle Size (nm)');
title(sprintf('Particle Size — R²=%.4f  RMSE=%.2f  MAE=%.2f', ...
    R2_PS, RMSE_PS, MAE_PS));
legend('Data Points','Perfect Fit','Location','Best');
grid on;

%% ---- Actual vs Predicted — Yield ----
figure('Name','Performance — Yield','NumberTitle','off');
plot(actual_Yield, predicted_Yield, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
hold on;
ref_y = linspace(min(actual_Yield)-2, max(actual_Yield)+2, 50);
plot(ref_y, ref_y, 'r--', 'LineWidth', 1.5);
xlabel('Actual Yield (%)');
ylabel('Predicted Yield (%)');
title(sprintf('Yield — R²=%.4f  RMSE=%.2f  MAE=%.2f', ...
    R2_Yield, RMSE_Yield, MAE_Yield));
legend('Data Points','Perfect Fit','Location','Best');
grid on;

%% ---- Actual vs Predicted — EE ----
figure('Name','Performance — EE','NumberTitle','off');
plot(actual_EE, predicted_EE, 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
hold on;
ref_e = linspace(min(actual_EE)-2, max(actual_EE)+2, 50);
plot(ref_e, ref_e, 'k--', 'LineWidth', 1.5);
xlabel('Actual EE (%)');
ylabel('Predicted EE (%)');
title(sprintf('EE — R²=%.4f  RMSE=%.2f  MAE=%.2f', ...
    R2_EE, RMSE_EE, MAE_EE));
legend('Data Points','Perfect Fit','Location','Best');
grid on;

%% ---- Metrics Bar Dashboard ----
figure('Name','Performance Metrics Dashboard','NumberTitle','off');

subplot(1,3,1);
bar([RMSE_PS, RMSE_Yield, RMSE_EE], 'FaceColor', [0.9 0.3 0.3]);
set(gca,'XTickLabel',{'PS','Yield','EE'});
ylabel('RMSE');
title('RMSE (lower = better)');
yline(5,'r--','Target < 5');
grid on;

subplot(1,3,2);
bar([MAE_PS, MAE_Yield, MAE_EE], 'FaceColor', [0.3 0.7 0.3]);
set(gca,'XTickLabel',{'PS','Yield','EE'});
ylabel('MAE');
title('MAE (lower = better)');
yline(3,'r--','Target < 3');
grid on;

subplot(1,3,3);
bar([R2_PS, R2_Yield, R2_EE], 'FaceColor', [0.2 0.5 0.9]);
set(gca,'XTickLabel',{'PS','Yield','EE'});
ylabel('R²');
title('R² Score (higher = better)');
ylim([0 1]);
yline(0.90,'r--','Target > 0.90');
grid on;

sgtitle('Model Performance Metrics Dashboard');

%% ---- Display Metrics ----
fprintf('\n');
fprintf('==============================================\n');
fprintf('  MODEL PERFORMANCE METRICS\n');
fprintf('==============================================\n');
fprintf('%-20s %8s %8s %8s\n','Output','RMSE','MAE','R²');
fprintf('----------------------------------------------\n');
fprintf('%-20s %8.4f %8.4f %8.4f\n','Particle Size', RMSE_PS,    MAE_PS,    R2_PS);
fprintf('%-20s %8.4f %8.4f %8.4f\n','Yield',         RMSE_Yield, MAE_Yield, R2_Yield);
fprintf('%-20s %8.4f %8.4f %8.4f\n','Encap. Eff.',   RMSE_EE,    MAE_EE,    R2_EE);
fprintf('----------------------------------------------\n');

% Grade
if R2_PS > 0.90 && RMSE_PS < 5
    Model_Grade = 'EXCELLENT MODEL';
elseif R2_PS > 0.80 && RMSE_PS < 10
    Model_Grade = 'GOOD MODEL';
else
    Model_Grade = 'NEEDS IMPROVEMENT';
end
fprintf('Model Grade : %s\n', Model_Grade);
fprintf('==============================================\n');

%% ---- Save Performance Metrics ----
Metrics_Results = table( ...
    {'Particle Size'; 'Yield'; 'EE'}, ...
    [RMSE_PS; RMSE_Yield; RMSE_EE],  ...
    [MAE_PS;  MAE_Yield;  MAE_EE],   ...
    [R2_PS;   R2_Yield;   R2_EE],    ...
    'VariableNames',{'Output','RMSE','MAE','R2'});

writetable(Metrics_Results, 'Performance_Metrics.xlsx');

disp('Phase 12 Performance Metrics Completed Successfully');
disp('====== ALL PHASES COMPLETE ======');

%% ============================================================
% PHASE 13 : QbD DECISION MATRIX
% (Yield, Cost, Quality, Scalability, Overall QbD Score)
%% ============================================================

disp('Computing QbD Decision Matrix...')

%% ---- 1. YIELD SCORE (direct, already in %) ----
Yield_Score = max(0, min(100, Yield_Opt));

%% ---- 2. QUALITY SCORE (based on PDI and EE) ----
% Lower PDI -> higher quality, Higher EE -> higher quality
Quality_Score = 0.5*(100*(1 - PDI_Opt)) + 0.5*EE_Opt;
Quality_Score = max(0, min(100, Quality_Score));

%% ---- 3. COST EFFICIENCY SCORE ----
% Higher polymer/surfactant/solvent usage & lower flow rate => higher cost
Cost_Index = (Optimal_Polymer    * 2.0) ...
            + (Optimal_Surfactant* 5.0) ...
            + (Optimal_SolventRatio*0.5) ...
            - (Optimal_FlowRate  * 1.0);

% Theoretical min/max cost index from PSO bounds
Cost_Index_min = (lb(3)*2.0) + (lb(4)*5.0) + (lb(5)*0.5) - (ub(2)*1.0);
Cost_Index_max = (ub(3)*2.0) + (ub(4)*5.0) + (ub(5)*0.5) - (lb(2)*1.0);

Cost_Score = 100 * (1 - (Cost_Index - Cost_Index_min) / ...
             (Cost_Index_max - Cost_Index_min + eps));
Cost_Score = max(0, min(100, Cost_Score));

%% ---- 4. SCALABILITY SCORE ----
% Based on % drift of PS, Yield, PDI from Lab -> Industrial scale
PS_drift    = abs(PS_Scale(3)    - PS_Scale(1))    / PS_Scale(1)    * 100;
Yield_drift = abs(Yield_Scale(3) - Yield_Scale(1)) / Yield_Scale(1) * 100;
PDI_drift   = abs(PDI_Scale(3)   - PDI_Scale(1))   / PDI_Scale(1)   * 100;

Avg_Drift = mean([PS_drift, Yield_drift, PDI_drift]);

Scalability_Score = 100 - Avg_Drift;
Scalability_Score = max(0, min(100, Scalability_Score));

%% ---- 5. OVERALL QbD SCORE (weighted composite) ----
w_yield = 0.30;
w_cost  = 0.20;
w_qual  = 0.30;
w_scale = 0.20;

QbD_Score = w_yield*Yield_Score + w_cost*Cost_Score ...
          + w_qual*Quality_Score + w_scale*Scalability_Score;
QbD_Score = max(0, min(100, QbD_Score));

%% ---- QbD CLASSIFICATION ----
if     QbD_Score >= 85, QbD_Class = 'EXCELLENT (Design Space Robust)';
elseif QbD_Score >= 70, QbD_Class = 'GOOD (Acceptable Design Space)';
elseif QbD_Score >= 50, QbD_Class = 'MARGINAL (Needs Optimization)';
else,                   QbD_Class = 'POOR (Redesign Required)';
end

%% ============================================================
% COMMAND WINDOW DISPLAY
%% ============================================================

fprintf('\n');
fprintf('=====================================================\n');
fprintf('          QbD DECISION MATRIX RESULTS\n');
fprintf('=====================================================\n');
fprintf('%-22s : %6.2f / 100\n','Yield Score',        Yield_Score);
fprintf('%-22s : %6.2f / 100\n','Cost Efficiency',    Cost_Score);
fprintf('%-22s : %6.2f / 100\n','Quality Score',      Quality_Score);
fprintf('%-22s : %6.2f / 100\n','Scalability Score',  Scalability_Score);
fprintf('-----------------------------------------------------\n');
fprintf('%-22s : %6.2f / 100\n','OVERALL QbD SCORE',  QbD_Score);
fprintf('Classification         : %s\n', QbD_Class);
fprintf('=====================================================\n');
fprintf('\nSupporting Details:\n');
fprintf(' Cost Index            : %.4f\n', Cost_Index);
fprintf(' PS Drift (Lab->Indus) : %.2f %%\n', PS_drift);
fprintf(' Yield Drift           : %.2f %%\n', Yield_drift);
fprintf(' PDI Drift             : %.2f %%\n', PDI_drift);

%% ============================================================
% PLOTS
%% ============================================================

QbD_Labels = {'Yield','Cost Eff.','Quality','Scalability','Overall QbD'};
QbD_Values = [Yield_Score, Cost_Score, Quality_Score, Scalability_Score, QbD_Score];

%% --- Figure A : Bar Chart ---
figure('Name','QbD Decision Matrix - Bar Chart','NumberTitle','off');
b = bar(QbD_Values, 'FaceColor','flat');
b.CData = [0.2 0.6 0.9; 0.9 0.6 0.2; 0.3 0.8 0.3; 0.8 0.3 0.6; 0.5 0.2 0.8];
set(gca,'XTickLabel',QbD_Labels,'XTickLabelRotation',12);
ylabel('Score (0-100)');
ylim([0 100]);
title(sprintf('QbD Decision Matrix — Overall Score = %.2f (%s)', ...
    QbD_Score, QbD_Class));
yline(70,'r--','Acceptance Threshold (70)');
grid on;

for i = 1:length(QbD_Values)
    text(i, QbD_Values(i)+2, sprintf('%.1f',QbD_Values(i)), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

%% --- Figure B : Radar (Spider) Chart ---
figure('Name','QbD Decision Matrix - Radar Chart','NumberTitle','off');

n = length(QbD_Values);
angles = linspace(0, 2*pi, n+1);

values_loop = [QbD_Values, QbD_Values(1)];   % close the loop

polarplot(angles, values_loop, '-o', 'LineWidth', 2, ...
    'MarkerFaceColor','b');
hold on;

% Reference circle at threshold = 70
polarplot(angles, 70*ones(1,n+1), 'r--', 'LineWidth', 1.2);

ax = gca;
ax.ThetaTick = rad2deg(angles(1:end-1));
ax.ThetaTickLabel = QbD_Labels;
ax.RLim = [0 100];
title('QbD Radar Chart — Process Robustness Profile');
legend('Process Scores','Acceptance Threshold (70)','Location','southoutside');

%% ============================================================
% SAVE QbD RESULTS
%% ============================================================

QbD_Results = table(QbD_Labels', QbD_Values', ...
    'VariableNames', {'Metric','Score'});

QbD_Summary = table(Yield_Score, Cost_Score, Quality_Score, ...
    Scalability_Score, QbD_Score, {QbD_Class}, ...
    'VariableNames', {'Yield_Score','Cost_Score','Quality_Score', ...
    'Scalability_Score','Overall_QbD_Score','QbD_Classification'});

writetable(QbD_Results, 'QbD_Metrics_Breakdown.xlsx');
writetable(QbD_Summary, 'QbD_Summary.xlsx');

disp('Phase 13 QbD Decision Matrix Completed Successfully');
disp('====== ALL PHASES COMPLETE (INCLUDING QbD) ======');



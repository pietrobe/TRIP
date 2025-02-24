clear;
restoredefaultpath;

% Path to the data folder
%data_path = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD/';
%data_path = "/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD.128_RII_CONTRIB_FAST/"
%addpath(data_path)

%% Strange behaviour
% profiles_CRD_7_12;

i_theta_v = [5 6 7 8];
%i_theta_v = [5];

X = 10;
Y = 40;
CRD_PRD = "PRD";

i_chi = 1;

path_prd_v = ["/work/sriva/JobsMareNostrumTRIP/out_Results/output/main/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_fix_conv_KQ_MC.pmd.PRD" ]
     % "/home/sriva/ms5prj/output/test_RTOL/1e-7/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd.PRD" ];

% path_crd = "/home/sriva/hero_scratch/dati3D/testTRIP/cai_1Bx_1By_1Bz_1Vx_1Vy_1Vz_GT4_32x32x133.pmd.CRD";
% path_crd_2 = "/home/sriva/hero_scratch/dati3D/testTRIP/orig/output_rewind_code/";

% path_crd = "/work/sriva/JobsMareNostrumTRIP/out_Results/output/main/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd.CRD_limit";
% path_crd = "/work/sriva/mu1/AR_385_Cut_64x64_mirrorxy-CRD_I_V0-B0_V0_conv_KQ.pmd.CRD/";
% path_crd_2 = "";

%path_n = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD';
%path_f = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd_RII_CONTRIB_FAST';


M = "none"; % Marker
M_vec = ["x" "o"];
M_Size = [8 12];

% Use one of these values: '+' | 'o' | '*' | '.' | 'x' |
% 'square' | 'diamond' | 'v' | '^' | '>' | '<' |
% 'pentagram' | 'hexagram' | '|' |
% '_' | 'none'.

for ii = 1:size(path_prd_v,1)
    
    path_prd = path_prd_v(ii);
    fprintf("in path (%d): %s\n", ii, path_prd);

    for i_theta = i_theta_v
        plot_profiles(path_prd, ...
            i_theta, i_chi, ...
            X, Y, CRD_PRD, ...
            M_vec(ii), M_Size(ii));

        % plot_profiles(path_prd_2, ...
        %     i_theta, i_chi, "o", 10);

    end
end

return;

% Add Porta result

file_porta = "/work/sriva/JobsMareNostrumTRIP/out_Results/resharpintensity/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_KQ_MC_quad1_(10-10).txt";
add_Porta_result(file_porta);

%% Functions

function add_Porta_result(file_name)

t = readmatrix("/work/sriva/JobsMareNostrumTRIP/out_Results/resharpintensity/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_KQ_MC_quad1_(10-10).txt");
fp =  vacuum_to_air(t(:,1));
subplot(2,2,1)
plot(fp , t(:,2), Marker="o")
hold on

subplot(2,2,2)
plot(fp , -100*t(:,3)./t(:,2), Marker="o")
hold on

subplot(2,2,3)
plot(fp , -100*t(:,4)./t(:,2), Marker="o")
hold on

subplot(2,2,4)
plot(fp , 100*t(:,5)./t(:,2), Marker="o", DisplayName="Porta")

hold on

legend

end

%% plot profiles
function plot_profiles(data_path, i_theta, i_chi, ...
    X, Y, ...
    CRD_PRD, ...
    M, MS)

addpath(data_path)
cmd_profiles = sprintf("profiles_%s_%d_%d", CRD_PRD,  X, Y);

fprintf("in: %s\n", cmd_profiles);

run(cmd_profiles)

% profiles_PRD_mu0960290_2_2;

if i_theta <= size(Field,2)/2
    disp('WARNING: mu < 0! Setting i_theta to the first emerging direction:')

    i_theta = size(Field,2)/2 + 1;

    disp(i_theta)
end

% load_nu_grid;
% nu_grid = nu_grid_(:,1);
nu_grid = vacuum_to_air(freq_to_angstrom(nu_grid_));
% nu_grid = freq_to_angstrom(nu_grid_);

mu  = mu_grid(i_theta);
chi = chi_grid(i_chi);

% disp(mu_grid);

fprintf("mu = %f, acos nu = %f = %f deg \n", mu, acos(mu), rad2deg(acos(mu)));
fprintf("chi = %f = %f deg\n", chi, rad2deg(chi));

xmin = 4220.75;
xmax = 4233.1;

update_plot(nu_grid, Field, i_theta, i_chi, mu, chi, xmin, xmax, M, MS)

sgtitle(sprintf("Coordinates X: %d, Y: %d", X, Y))

rmpath(data_path)
end


%% Update plot
function update_plot(nu_grid, Field, i_theta, i_chi, mu, chi, xmin, xmax, M, MS)
subplot(2,2,1)
set_plot_defaults
plot(nu_grid,Field{1,i_theta,i_chi}, Marker=M, MarkerSize=MS)
xlim([xmin, xmax])
hold on
title(['$(\mu,\chi) = ($', num2str(mu,3), ', ', num2str(chi,3),')'])
ylabel('$I$')
set_plot_defaults

subplot(2,2,2)
set_plot_defaults
plot(nu_grid,Field{2,i_theta,i_chi}, Marker=M, MarkerSize=MS)
xlim([xmin, xmax])
hold on
ylabel('$Q/I[\%]$')
set_plot_defaults

subplot(2,2,3)
set_plot_defaults
plot(nu_grid,Field{3,i_theta,i_chi}, Marker=M, MarkerSize=MS)
xlim([xmin, xmax])
hold on
ylabel('$U/I[\%]$')
set_plot_defaults

subplot(2,2,4)
%DisplayName=sprintf("mu=%lf", mu), ...
set_plot_defaults
plot(nu_grid,Field{4,i_theta,i_chi}, ...
    'DisplayName', sprintf("$\\mu=%.2f$", mu), ...
    Marker=M, MarkerSize=MS);

legend;
xlim([xmin, xmax])
hold on
ylabel('$V/I[\%]$')
set_plot_defaults


end

function set_plot_defaults()

set(gca,'FontSize', 16);
set(gca,'TickLabelInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 2);
grid on

end



% Function to convert centimeters to angstroms
function angstrom = cm_to_angstrom(cm)
angstrom = cm .* 1e8;
end


% Function to convert frequency to angstroms
function angstrom = freq_to_angstrom(freq)
speed_of_light = 29979245800; % Speed of light in centmeters per second
cm = speed_of_light ./ freq;
angstrom = cm_to_angstrom(cm);
end

% Function to convert vacuum wavelength to air wavelength
function wave_air = vacuum_to_air(wave, to_air_limit)
% Default value for to_air_limit if not provided
if nargin < 2
    to_air_limit = 200.0;
end

wave2 = (wave .* wave) / 100.0;

fact = 1.0 + 2.735182e-4 + (1.314182e0 + (1.0 / 2.76249e+4) * wave2) ./ wave2;
fact = fact .* (wave > to_air_limit) + 1.0 * (wave < to_air_limit);

wave_air = wave ./ fact;
end

function cont_map = make_continuum_map(model_size, i_theta, i_chi, input_dir, file_base)

cont_map = zeros(model_size, model_size);

for i = 1:model_size
    for j = 1:model_size
        input_script = sprintf(file_base, (i-1), (j-1));
        run( fullfile(input_dir, input_script));

        field_I = Field{1,i_theta,i_chi};
        cont_map(i, j) = field_I(1);
    end

    if i == 1
        mu  = mu_grid(i_theta);
        chi = chi_grid(i_chi);

        fprintf("mu = %f, acos nu = %f = %f deg \n", mu, acos(mu), rad2deg(acos(mu)));
        fprintf("chi = %f = %f deg\n", chi, rad2deg(chi));
    end

    if mod(i, 4) == 0
        fprintf("i = %d \n", i);
    end
end
end

function cont_map = make_continuum_map_nu1(model_size, i_theta, i_chi, input_dir, file_base)

cont_map = zeros(model_size, model_size);

for i = 1:model_size
    for j = 1:model_size
        input_script = sprintf(file_base, (i-1), (j-1));
        run( fullfile(input_dir, input_script));

        field_I = Field_Omega{1};
        cont_map(i, j) = field_I(1);
    end

    % if i == 1
    %     mu  = mu_grid(i_theta);
    %     chi = chi_grid(i_chi);
    %
    %     fprintf("mu = %f, acos nu = %f = %f deg \n", mu, acos(mu), rad2deg(acos(mu)));
    %     fprintf("chi = %f = %f deg\n", chi, rad2deg(chi));
    % end

    if mod(i, 4) == 0
        fprintf("i = %d \n", i);
    end
end
end





%
% plot_profiles(path_crd, ...
%     i_theta, i_chi, "*");

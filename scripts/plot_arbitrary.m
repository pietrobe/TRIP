clear all;
% close all;


input_path = "/home/sriva/ms5prj/output/test_RTOL/1e-9/64x64/AR_385_Cut_64x64_mirrorxy-CRD_I_V0_fix_conv_KQ_MC.pmd.PRD";

mu_v = ["01000" "01834" "03000" "07000" "09603" "09900" "10000"];
mu_v = ["01000" "03000" "07000" "09900" "10000"];

% mu_v = [ "01834" ];

chi = "01963";
chi = "00000";

X = 42;
Y = 42;

PRD_CRD = "PRD";

for mu = mu_v
    plot_arbitrary_mu(input_path, PRD_CRD, mu, chi, X, Y, ...
        "O", 3)
end

%% function plot_arbitrary_mu
function plot_arbitrary_mu(input_dir, PRD_CRD, ...
    mu, chi, ...
    X, Y, ...
    M, MSize)

mu_ = str2double(mu) / 10000.0;
chi_ = str2double(chi) / 10000.0;

input_file = sprintf("profiles_%s_mu%s_chi%s_%d_%d", PRD_CRD, mu, chi, X, Y);
input_file_b = sprintf("profiles_%s_%d_%d", PRD_CRD, X, Y);

full_infile = fullfile(input_dir, input_file);
full_infile_b = fullfile(input_dir, input_file_b);

run(full_infile);
run(full_infile_b);
nu_grid = vacuum_to_air(freq_to_angstrom(nu_grid_));

subplot(2,2,1)
plot(nu_grid, Field_Omega{1}, Marker=M, MarkerSize=MSize);
title("I");
grid on
hold on
set_plot_defaults();

subplot(2,2,2)
plot(nu_grid, Field_Omega{2}, Marker=M, MarkerSize=MSize);
grid on
title("Q/I");
hold on
set_plot_defaults();

subplot(2,2,3)
plot(nu_grid, Field_Omega{3}, Marker=M, MarkerSize=MSize);
grid on
title("U/I");
hold on
set_plot_defaults();

subplot(2,2,4)
plot(nu_grid, Field_Omega{4}, ...
    DisplayName=sprintf("$\\mu=%.3f$", mu_), ...
    Marker=M, MarkerSize=MSize);
% legend(sprintf("$\\mu=%s$", mu_), 'Interpreter', 'latex');
legend
grid on
title("V/I");
hold on

set_plot_defaults();

sgtitle(sprintf("Coordinates X: %d, Y: %d", X, Y))

end

%% set_plot_defaults
function set_plot_defaults()

set(gca,'FontSize', 16);
set(gca,'TickLabelInterpreter', 'Latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLineLineWidth', 2);
grid on

end

%% Function to convert vacuum wavelength to air wavelength
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
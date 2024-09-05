clear all;

% Path to the data folder
data_path = '/home/sriva/hero_scratch/Comparison-TRIP-PORTA/PORTA/output/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd/';  

addpath(data_path)

profiles_CRD_1_1

i_theta = 8;
i_chi = 1;

if i_theta <= size(Field,2)/2
    disp('WARNING: mu < 0! Setting i_theta to the first emerging direction:')    

    i_theta = size(Field,2)/2 + 1;

    disp(i_theta)

end

% load_nu_grid;
% nu_grid = nu_grid_(:,1);
nu_grid = vacuum_to_air(freq_to_angstrom(nu_grid_));

mu  = mu_grid(i_theta);
chi = chi_grid(i_chi);

xmin = 4225;
xmax = 4228;

M = "none"; % Marker
% Use one of these values: '+' | 'o' | '*' | '.' | 'x' |
% 'square' | 'diamond' | 'v' | '^' | '>' | '<' | 
% 'pentagram' | 'hexagram' | '|' |
% '_' | 'none'.


subplot(2,2,1)
set_plot_defaults
plot(nu_grid,Field{1,i_theta,i_chi}, Marker=M)
xlim([xmin, xmax])
hold on
title(['$(\mu,\chi) = ($', num2str(mu,3), ', ', num2str(chi,3),')'])
ylabel('$I$')
set_plot_defaults

subplot(2,2,2)
set_plot_defaults
plot(nu_grid,Field{2,i_theta,i_chi}, Marker=M)
xlim([xmin, xmax])
hold on
ylabel('$Q/I[\%]$')
set_plot_defaults

subplot(2,2,3)
set_plot_defaults
plot(nu_grid,Field{3,i_theta,i_chi}, Marker=M)
xlim([xmin, xmax])
hold on
ylabel('$U/I[\%]$')
set_plot_defaults

subplot(2,2,4)
set_plot_defaults
plot(nu_grid,Field{4,i_theta,i_chi}, Marker=M)
xlim([xmin, xmax])
hold on
ylabel('$V/I[\%]$')
set_plot_defaults

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



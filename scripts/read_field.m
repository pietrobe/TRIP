clear;
restoredefaultpath;

% Path to the data folder
%data_path = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD/';
%data_path = "/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD.128_RII_CONTRIB_FAST/"
%addpath(data_path)

i_theta = 8;
i_chi = 1;

%path_prd = "/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD.1000G.RII";
path_prd = "/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD.1000G.RII_FAST/"

%path_crd = "/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.CRD_test";

%path_n = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd.PRD';
%path_f = '/home/sriva/hero_scratch/TRIP_test/cai_0Bx_0By_0Bz_0Vx_0Vy_0Vz_GT4_5x5x133_it100.pmd_RII_CONTRIB_FAST';


plot_profiles(path_prd, ...
    i_theta, i_chi)

function plot_profiles(data_path, i_theta, i_chi)
    addpath(data_path)
    
    profiles_PRD_4_4;
    
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
    
    xmin = 4224.5;
    xmax = 4228.5;
    
    M = "o"; % Marker
    % Use one of these values: '+' | 'o' | '*' | '.' | 'x' |
    % 'square' | 'diamond' | 'v' | '^' | '>' | '<' | 
    % 'pentagram' | 'hexagram' | '|' |
    % '_' | 'none'.
    
    update_plot(nu_grid, Field, i_theta, i_chi, mu, chi, xmin, xmax, M)

    rmpath(data_path)
end



function update_plot(nu_grid, Field, i_theta, i_chi, mu, chi, xmin, xmax, M)
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



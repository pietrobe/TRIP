function [] = plot_ar(results, varargin)
%PLOT_AR Plot emergent field results
%
% plot_ar(results)
% plot_ar(results, 'param', value, ...)
%
% Inputs:
%   results - struct or array of structs from read_ar_emergent_field
%             Each struct should have: I, QI, UI, VI, frequencies, mu_used, chi_used
%
% Options:
%   'title'    - plot title (default: 'Emergent Field')
%   'wavelength' - if true, convert freq to wavelength in air (default: true)

p = inputParser;
addParameter(p, 'title', 'Emergent Field');
addParameter(p, 'wavelength', true);
parse(p, varargin{:});

if ~iscell(results)
    results = {results};
end

n_results = length(results);

figure;
for idx = 1:n_results
    r = results{idx};
    
    if p.Results.wavelength
        wave_vac = freq_to_angstrom(r.frequencies);
        x = vacuum_to_air(wave_vac);
        xlabel_str = 'Wavelength (Angstrom, air)';
    else
        x = r.frequencies;
        xlabel_str = 'Frequency (Hz)';
    end
    
    QI_pct = (r.QI);
    UI_pct = (r.UI);
    VI_pct = (r.VI);
    
    subplot(2, 2, 1);
    plot(x, r.I, '-');
    hold on;
    xlabel(xlabel_str);
    ylabel('Stokes I');
    title('Stokes I');
    if n_results > 1
        legend_str = cell(1, n_results);
        for k = 1:n_results
            legend_str{k} = sprintf('\\mu=%.3f, \\chi=%.3f', results{k}.mu_used, results{k}.chi_used);
        end
        legend(legend_str, 'Location', 'best');
    else
        legend(sprintf('\\mu=%.3f, \\chi=%.3f', r.mu_used, r.chi_used), 'Location', 'best');
    end
    grid on;
    
    subplot(2, 2, 2);
    plot(x, QI_pct, '-');
    hold on;
    xlabel(xlabel_str);
    ylabel('Q/I (%)');
    title('Stokes Q/I');
    grid on;
    
    subplot(2, 2, 3);
    plot(x, UI_pct, '-');
    hold on;
    xlabel(xlabel_str);
    ylabel('U/I (%)');
    title('Stokes U/I');
    grid on;
    
    subplot(2, 2, 4);
    plot(x, VI_pct, '-');
    hold on;
    xlabel(xlabel_str);
    ylabel('V/I (%)');
    title('Stokes V/I');
    grid on;
end

sgtitle(p.Results.title);
end
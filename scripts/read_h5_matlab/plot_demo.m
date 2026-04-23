function [] = plot_demo(i, j, incl, az, file_name, file_ref)
%PLOT_DEMO Summary of this function goes here
%   Detailed explanation goes here

% Add path to utility functions
script_path = '.';
addpath(script_path);

% time the HDF5 read only
tRead = tic;
[nu_grid__, mu_grid_, theta_grid_, chi_grid_, Field_] = THDF_read_surface_field_from_hdf5(file_name, i, j);
readElapsed = toc(tRead);
fprintf('THDF_read_surface_field_from_hdf5: elapsed = %.3f s\n', readElapsed);

% Convert frequency (Hz) to wavelength in air (Angstrom)
nu_grid_vacuum_angstrom = freq_to_angstrom(nu_grid__);
wave_air = vacuum_to_air(nu_grid_vacuum_angstrom);
wave_air = wave_air(:);

fprintf('wave_air (1:10): '); disp(wave_air(1:10).');

% execute the reference script to get reference data if provided
% treat missing, empty, or empty-string `file_ref` as not provided
ref_present = (nargin >= 6) && ~isempty(file_ref) && ~(isstring(file_ref) && (strlength(file_ref) == 0 || strcmpi(file_ref, 'none')));
if ref_present
    refRead = tic;
	run(char(file_ref));
    refElapsed = toc(refRead);
    fprintf('Reference script execution: elapsed = %.3f s\n', refElapsed);
end

if ref_present
	fprintf('Reference data loaded from %s\n', char(file_ref));
else
	fprintf('No reference data provided.\n');
end

figure;
subplot(2,2,1);
dataI = Field_{1, incl, az};
if ref_present
	wave_air_ref = vacuum_to_air(freq_to_angstrom(nu_grid_));
	plot(wave_air, dataI, 'b-', wave_air_ref, Field{1,incl,az}, 'ro');
	legend('Read from HDF5', 'Reference data'); 
else
	plot(wave_air, dataI, 'b-');
	legend('Read from HDF5');
end
xlabel('Wavelength (Angstrom, air)');
ylabel('Stokes I');
title('Stokes I comparison');

subplot(2,2,2);
if ref_present
	plot(wave_air, Field_{2, incl, az}, 'b-', wave_air_ref, Field{2,incl,az}, 'ro');
	legend('Read from HDF5', 'Reference data');
else
	plot(wave_air, Field_{2, incl, az}, 'b-');
	legend('Read from HDF5');
end
xlabel('Wavelength (Angstrom, air)');
ylabel('Stokes Q');
title('Stokes Q comparison');

subplot(2,2,3);
if ref_present
	plot(wave_air, Field_{3, incl, az}, 'b-', wave_air_ref, Field{3,incl,az}, 'ro');
	legend('Read from HDF5', 'Reference data');
else
	plot(wave_air, Field_{3, incl, az}, 'b-');
	legend('Read from HDF5');
end
xlabel('Wavelength (Angstrom, air)');
ylabel('Stokes U');
title('Stokes U comparison');

subplot(2,2,4);
if ref_present
	plot(wave_air, Field_{4, incl, az}, 'b-', wave_air_ref, Field{4,incl,az}, 'ro');
	legend('Read from HDF5', 'Reference data');
else
	plot(wave_air, Field_{4, incl, az}, 'b-');
	legend('Read from HDF5');
end
xlabel('Wavelength (Angstrom, air)');
ylabel('Stokes V');
title('Stokes V comparison');
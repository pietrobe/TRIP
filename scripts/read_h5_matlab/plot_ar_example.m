clear all;
close all;

addpath('/home/simone/git/hdf5_demo/output/matlab');

file_name_ar = '/home/simone/cscs_scratch/PORTA_OUT/h5_Borig_test88_1/small_ssd_from15D_8x8x143.h5.PRD/emergent_field_abitrary_Omega.h5';

i = 4;
j = 4;

fprintf('Reading mu=1.0, chi=0.0...\n');
result1 = read_ar_emergent_field(file_name_ar, i, j, 0.1, 0.0);

fprintf('\nReading mu=0.3, chi=0.0...\n');
result2 = read_ar_emergent_field(file_name_ar, i, j, 0.3, 0.0);

fprintf('\nPlotting...\n');
plot_ar({result1, result2}, 'title', 'Emergent Field Comparison');
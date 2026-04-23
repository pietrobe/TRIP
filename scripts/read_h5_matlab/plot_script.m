clear all;
close all;

i = 4;
j = 4;

incl = 1;
az = 1;

base_out_path = "/run/host/home/simone/cscs_scratch/PORTA_OUT/h5_Borig_test88_1/small_ssd_from15D_8x8x143.h5.PRD";

% file_ref = ""; % disable reference for now

file_name = fullfile(base_out_path, "emergent_field_angular_grid.h5");
file_name_ar = fullfile(base_out_path, "emergent_field_abitrary_Omega.h5");

file_ref = sprintf(file_name, i-1, j-1);
fprintf("%s\n", file_name);

% file_ref = ""; % disable reference for now

% file_ref = ""; % disable reference for now

plot_demo(i, j, incl, az, file_name, "");

result = read_ar_emergent_field(file_name_ar, i, j, 1, 0);
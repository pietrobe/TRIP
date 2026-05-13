function result = read_ar_emergent_field(h5file, i, j, mu_req, chi_req)
%READ_AR_EMERGENT_FIELD Read emergent field at given position and angle
% Read the emergent Stokes parameters (I, Q/I, U/I, V/I) for a specific spatial position (i,j) and beam direction (mu, chi) from the HDF5 file. 
% The function finds the closest available beam direction to the requested (mu_req, chi_req) 
% and reads the corresponding data.
% result = read_ar_emergent_field(h5file, i, j, mu_req, chi_req)
%
% Inputs:
%   h5file   - path to HDF5 file
%   i, j     - spatial indices (1-based)
%   mu_req   - requested mu value (cos of inclination angle)
%   chi_req  - requested chi value (azimuthal angle)
%
% Outputs:
%   result - struct with fields:
%     .I  - Stokes I (column vector, 96 frequencies)
%     .QI - Stokes Q (column vector)
%     .UI - Stokes U (column vector)
%     .VI - Stokes V (column vector)
%     .frequencies - frequency vector (Hz)
%     .mu_used - mu value actually used
%     .chi_used - chi value actually used

if nargin < 5
    error('Usage: read_ar_emergent_field(h5file, i, j, mu_req, chi_req)');
end

h5file = char(h5file);

mu_all = h5read(h5file, '/beams_directions/mu');
chi_all = h5read(h5file, '/beams_directions/chi');
N_dirs = h5read(h5file, '/beams_directions/N_directions');

mu_all = double(mu_all(:));
chi_all = double(chi_all(:));
N_dirs = double(N_dirs);

dist = sqrt((mu_all - mu_req).^2 + (chi_all - chi_req).^2);
[~, idx] = min(dist);

mu_used = mu_all(idx);
chi_used = chi_all(idx);

if i == 0, i = 1; fprintf('Converted i=0 to 1\n'); end
if j == 0, j = 1; fprintf('Converted j=0 to 1\n'); end

fprintf('Requested mu=%.3f, chi=%.3f\n', mu_req, chi_req);
fprintf('Using mu=%.3f, chi=%.3f (dir index %d)\n', mu_used, chi_used, idx);

ds_name = sprintf('mu_%.3f_chi_%.3f', mu_used, chi_used);

ds_path_I = sprintf('/emergent_beams/beams_I/%s', ds_name);
ds_path_QI = sprintf('/emergent_beams/beams_QI_pc/%s', ds_name);
ds_path_UI = sprintf('/emergent_beams/beams_UI_pc/%s', ds_name);
ds_path_VI = sprintf('/emergent_beams/beams_VI_pc/%s', ds_name);

infoI = h5info(h5file, ds_path_I);
dims = double(infoI.Dataspace.Size);
fprintf('Dataset shape: [%s]\n', num2str(dims));

dims = dims(end:-1:1);
fprintf('Canonical shape (x,y,freq): [%s]\n', num2str(dims));

N_freq = dims(3);

start_idx = [i, j, 1];
count_idx = [1, 1, N_freq];

fprintf('Reading with start=[%s], count=[%s]\n', num2str(start_idx), num2str(count_idx));

try
    dataI = h5read(h5file, ds_path_I, start_idx, count_idx);
    fprintf('dataI after read: size=%s\n', num2str(size(dataI)));
    dataI = squeeze(dataI);
    fprintf('dataI after squeeze: size=%s\n', num2str(size(dataI)));
    dataQI = h5read(h5file, ds_path_QI, start_idx, count_idx);
    dataQI = squeeze(dataQI);
    dataUI = h5read(h5file, ds_path_UI, start_idx, count_idx);
    dataUI = squeeze(dataUI);
    dataVI = h5read(h5file, ds_path_VI, start_idx, count_idx);
    dataVI = squeeze(dataVI);
catch ME1
    fprintf('Direct read failed: %s\n', ME1.message);
    fprintf('Trying reversed order...\n');
    start_rev = [1, 1, 1];
    count_rev = [N_freq, dims(2), dims(1)];
    fprintf('Reversed read: start=[%s], count=[%s]\n', num2str(start_rev), num2str(count_rev));
    dataI = h5read(h5file, ds_path_I, start_rev, count_rev);
    fprintf('dataI raw shape: %s\n', num2str(size(dataI)));
    dataI = permute(dataI, [3, 2, 1]);
    fprintf('dataI after permute: %s\n', num2str(size(dataI)));
    dataI = squeeze(dataI(i, j, :));
    dataQI = h5read(h5file, ds_path_QI, start_rev, count_rev);
    dataQI = permute(dataQI, [3, 2, 1]);
    dataQI = squeeze(dataQI(i, j, :));
    dataUI = h5read(h5file, ds_path_UI, start_rev, count_rev);
    dataUI = permute(dataUI, [3, 2, 1]);
    dataUI = squeeze(dataUI(i, j, :));
    dataVI = h5read(h5file, ds_path_VI, start_rev, count_rev);
    dataVI = permute(dataVI, [3, 2, 1]);
    dataVI = squeeze(dataVI(i, j, :));
end

result = struct( ...
    'I', double(dataI(:)), ...
    'QI', double(dataQI(:)), ...
    'UI', double(dataUI(:)), ...
    'VI', double(dataVI(:)), ...
    'frequencies', double(h5read(h5file, '/frequencies_grid/frequencies')), ...
    'mu_used', mu_used, ...
    'chi_used', chi_used ...
    );

fprintf('Output: I size=%d, QI size=%d, frequencies=%d\n', numel(result.I), numel(result.QI), numel(result.frequencies));
end
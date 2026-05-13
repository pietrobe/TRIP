function [nu_grid_, mu_grid, theta_grid, chi_grid, Field] = THDF_read_surface_field_from_hdf5(h5file, x_i, y_i)

%THDF_read_surface_field_from_hdf5 Read surface Stokes field for a given i,j
%
% [nu_grid_, theta_grid, chi_grid, Field] = THDF_read_surface_field_from_hdf5(h5file, x_i, y_i)
%
% Returns:
%  - nu_grid_: vector of frequencies (N_frequencies x 1)
%  - theta_grid: vector of inclination angles (N_inclination_angles x 1)
%  - chi_grid: vector of azimuthal angles (N_azimuthal_angles x 1)
%  - Field: cell array sized {4, N_inclination_angles, N_azimuthal_angles}
%           Field{p,ti,ci} is a column vector (N_frequencies x 1) for Stokes p:
%             p=1 -> I, p=2 -> Q, p=3 -> U, p=4 -> V
%
% Example:
%  [nu,theta,chi,Field] = THDF_read_surface_field_from_hdf5('output_field_mpi.h5', 1, 1)

% Map of HDF5 group and dataset names for easy modification

% Auto-detect naming convention by checking which groups/datasets exist
infoRoot = h5info(h5file);
rootNames = {infoRoot.Groups.Name};
if isfield(infoRoot, 'Datasets') && ~isempty(infoRoot.Datasets)
    rootDsNames = {infoRoot.Datasets.Name};
else
    rootDsNames = {};
end

% Check for specific groups (more precise than contains)
hasOutputFieldGroup = any(strcmp(rootNames, '/output_field'));
hasEmergentFieldGroup = any(strcmp(rootNames, '/emergent_field'));

% Also check datasets for backward compatibility
hasOutputFieldData = any(contains(rootDsNames, 'output_field_data'));
hasEmergentFieldData = any(contains(rootDsNames, 'emergent_'));

if hasOutputFieldGroup || hasOutputFieldData
    original_naming = true;
elseif hasEmergentFieldGroup || hasEmergentFieldData
    original_naming = false;
else
    error('Cannot detect HDF5 naming convention: neither /output_field nor /emergent_field found in root');
end

if original_naming
    HDF5_NAMES = struct( ...
        'angular_grid_func', 'THDF_read_angular_grid_from_hdf5', ...
        'output_field_group', '/output_field', ...
        'output_field_I', '/output_field/output_field_data_I', ...
        'output_field_Q', '/output_field/output_field_data_Q', ...
        'output_field_U', '/output_field/output_field_data_U', ...
        'output_field_V', '/output_field/output_field_data_V', ...
        'output_frequencies_group', '/output_frequencies_grid', ...
        'frequencies', 'frequencies' ...
        );
else
    HDF5_NAMES = struct( ...
        'angular_grid_func', 'THDF_read_angular_grid_from_hdf5', ...
        'output_field_group', '/emergent_field', ...
        'output_field_I', '/emergent_field/emergent_I', ...
        'output_field_Q', '/emergent_field/emergent_QI_pc', ...
        'output_field_U', '/emergent_field/emergent_UI_pc', ...
        'output_field_V', '/emergent_field/emergent_VI_pc', ...
        'output_frequencies_group', '/frequencies_grid', ...
        'frequencies', 'frequencies' ...
        );
end

if nargin < 3
    error('Usage: THDF_read_surface_field_from_hdf5(h5file, x_i, y_i)');
end
if ~(ischar(h5file) || isstring(h5file))
    error('h5file must be a filename (string).');
end
h5file = char(h5file);


% Read angular grid (inclination/azimuth and optional frequencies)
ang = feval(HDF5_NAMES.angular_grid_func, h5file);
theta_grid = ang.inclination_angles(:);
chi_grid = ang.azimuthal_angles(:);
N_incl = double(ang.N_inclination_angles);
N_azim = double(ang.N_azimuthal_angles);
nu_grid_ = double(ang.frequencies(:));
mu_grid = cos(theta_grid);


fprintf('Reading surface field at position x_i=%d, y_i=%d\n', x_i, y_i);
fprintf('  Number of inclination angles: %d\n', N_incl);
fprintf('  Number of azimuthal angles: %d\n', N_azim);
fprintf('  Number of frequencies: %d\n', length(nu_grid_));

% Try to get frequencies from angular reader; otherwise infer from dataset
nu_grid_ = ang.frequencies(:);

% Get dataset info for output_field I to infer dimensions
dsPathI = HDF5_NAMES.output_field_I;
try
    infoI = h5info(h5file, dsPathI);
catch
    error(['Dataset ''' dsPathI ''' not found in file']);
end
    dims = infoI.Dataspace.Size; % expected [x y inclam azim freq]
    if numel(dims) < 5
        error('Unexpected dataset dimensionality for %s', dsPathI);
    end
    dims = double(dims);
    % h5info may return dimensions in reverse (Fortran/C order differences).
    % Flip to canonical [x y incl azim freq] ordering.
    dims = dims(end:-1:1);

    % Allow user to pass zero-based indices (convert to 1-based), otherwise verify
    if x_i == 0
        x_i = 1;
        warning('Converted zero-based x_i to 1 (MATLAB is 1-based)');
    end
    if y_i == 0
        y_i = 1;
        warning('Converted zero-based y_i to 1 (MATLAB is 1-based)');
    end

    % Check provided i,j against dataset dimensions
    if x_i < 1 || x_i > dims(1) || y_i < 1 || y_i > dims(2)
        error('Requested indices x_i=%d,y_i=%d exceed dataset dims [%d %d ...].', x_i, y_i, dims(1), dims(2));
    end

    Nf = double(dims(5));
if isempty(nu_grid_)
    % Try reading /output_frequencies_grid/frequencies
    try
        freqs = h5read(h5file, [HDF5_NAMES.output_frequencies_group '/' HDF5_NAMES.frequencies]);
        nu_grid_ = double(freqs(:));
    catch
        % fallback to indices 1..Nf
        nu_grid_ = (1:double(Nf)).';
    end
end

    % Prepare hyperslab to read all angles for given x_i,y_i
    % h5read uses 1-based indexing for start. Dataset axis order may vary
    % so detect positions of [x y incl azim freq] inside `dims`.
    % dims currently is in canonical order after earlier flip attempt, but
    % some HDF5 builds return reversed ordering; we match by sizes.
    dataset_rank = numel(dims);
    if dataset_rank < 5
        error('Unexpected dataset rank: %d', dataset_rank);
    end

    % Known sizes
    want = struct('N_incl', double(N_incl), 'N_azim', double(N_azim), 'Nf', double(Nf));
    assigned = false(1, dataset_rank);
    pos_incl = find(dims == want.N_incl & ~assigned, 1);
    if isempty(pos_incl)
        pos_incl = find(abs(dims - want.N_incl) == min(abs(dims - want.N_incl)), 1);
    end
    assigned(pos_incl) = true;
    pos_azim = find(dims == want.N_azim & ~assigned, 1);
    if isempty(pos_azim)
        pos_azim = find(abs(dims - want.N_azim) == min(abs(dims - want.N_azim)), 1);
        if assigned(pos_azim)
            pos_azim = find(~assigned, 1);
        end
    end
    assigned(pos_azim) = true;
    pos_nf = find(dims == want.Nf & ~assigned, 1);
    if isempty(pos_nf)
        pos_nf = find(abs(dims - want.Nf) == min(abs(dims - want.Nf)), 1);
        if assigned(pos_nf)
            pos_nf = find(~assigned, 1);
        end
    end
    assigned(pos_nf) = true;

    % Remaining two positions correspond to x and y (in some order)
    rempos = find(~assigned);
    if numel(rempos) ~= 2
        error('Unable to determine x/y axes from dataset dims: [%s]', num2str(dims));
    end
    % Assign x,y by checking which remaining dimension can contain the requested indices
    if dims(rempos(1)) >= x_i && dims(rempos(2)) >= y_i
        pos_x = rempos(1); pos_y = rempos(2);
    elseif dims(rempos(2)) >= x_i && dims(rempos(1)) >= y_i
        pos_x = rempos(2); pos_y = rempos(1);
    else
        % fallback: assign in order
        pos_x = rempos(1); pos_y = rempos(2);
    end

    % Build start and count vectors in dataset axis order
    start_ds = ones(1, dataset_rank);
    count_ds = ones(1, dataset_rank);
    start_ds(pos_x) = double(x_i); count_ds(pos_x) = 1;
    start_ds(pos_y) = double(y_i); count_ds(pos_y) = 1;
    start_ds(pos_incl) = 1; count_ds(pos_incl) = double(N_incl);
    start_ds(pos_azim) = 1; count_ds(pos_azim) = double(N_azim);
    start_ds(pos_nf) = 1; count_ds(pos_nf) = double(Nf);

    % Debug print of mapping
    fprintf('  dataset dims: [%s]', num2str(dims));
    fprintf('  axis mapping (pos_x,pos_y,pos_incl,pos_azim,pos_nf): %d %d %d %d %d\n', pos_x, pos_y, pos_incl, pos_azim, pos_nf);
    fprintf('  start: [%s]\n', num2str(start_ds));
    fprintf('  count: [%s]\n', num2str(count_ds));

    % Use computed start/count to read
    start = start_ds;
    count = count_ds;

% Read I,Q,U,V datasets
dsQPath = HDF5_NAMES.output_field_Q;
dsUPath = HDF5_NAMES.output_field_U;
dsVPath = HDF5_NAMES.output_field_V;
% Try reading using computed start/count. If axis-order mismatch persists,
% try reversing start/count (some MATLAB/HDF5 builds expect reversed order).
readError = [];
try
    dataI = h5read(h5file, dsPathI, start, count);
    dataQ = h5read(h5file, dsQPath, start, count);
    dataU = h5read(h5file, dsUPath, start, count);
    dataV = h5read(h5file, dsVPath, start, count);
catch ME1
    readError = ME1;
    fprintf('  initial h5read failed: %s\n', ME1.message);
    % Attempt reversed-order start/count
    start_rev = start(end:-1:1);
    count_rev = count(end:-1:1);
    fprintf('  trying reversed start/count: start=[%s], count=[%s]\n', num2str(start_rev), num2str(count_rev));
    try
        dataI = h5read(h5file, dsPathI, start_rev, count_rev);
        dataQ = h5read(h5file, dsQPath, start_rev, count_rev);
        dataU = h5read(h5file, dsUPath, start_rev, count_rev);
        dataV = h5read(h5file, dsVPath, start_rev, count_rev);
        % If successful, note that we used reversed mapping
        fprintf('  reversed-order read succeeded; transposing data to canonical order.\n');
        % We must permute dimensions back to canonical [x y incl azim freq]
        % Determine permutation that maps reversed dims to canonical dims
        perm = numel(dims):-1:1;
        dataI = permute(dataI, perm);
        dataQ = permute(dataQ, perm);
        dataU = permute(dataU, perm);
        dataV = permute(dataV, perm);
    catch ME2
        % Both attempts failed
        error('Error reading output_field datasets: initial: %s ; reversed: %s', readError.message, ME2.message);
    end
end



fprintf('Reading surface field at position II x_i=%d, y_i=%d\n', x_i, y_i);
fprintf('  Dataset dimensions: [%d %d %d %d %d]\n', dims);
fprintf('  Number of inclination angles: %d\n', N_incl);
fprintf('  Number of azimuthal angles: %d\n', N_azim);
fprintf('  Number of frequencies: %d\n', length(nu_grid_));

% dataI is expected shape [1 1 N_incl N_azim Nf]
% Build Field cell array: {4, N_incl, N_azim}
Field = cell(4, N_incl, N_azim);
for ti = 1:N_incl
    for ci = 1:N_azim
        % Extract frequency vector and ensure column shape
        vI = squeeze(dataI(1,1,ti,ci,:)); vI = double(vI(:));
        vQ = squeeze(dataQ(1,1,ti,ci,:)); vQ = double(vQ(:));
        vU = squeeze(dataU(1,1,ti,ci,:)); vU = double(vU(:));
        vV = squeeze(dataV(1,1,ti,ci,:)); vV = double(vV(:));
        Field{1,ti,ci} = vI;
        Field{2,ti,ci} = vQ;
        Field{3,ti,ci} = vU;
        Field{4,ti,ci} = vV;
    end
end

theta_grid = [pi-flip(theta_grid); theta_grid];
mu_grid = cos(theta_grid);




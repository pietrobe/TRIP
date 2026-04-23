function out = THDF_read_angular_grid_from_hdf5(h5file)

% Read /emergent_angular_grid from an HDF5 file and return a struct with:
% - N_directions
% - N_inclination_angles
% - N_azimuthal_angles
% - inclination_angles (double vector)
% - azimuthal_angles (double vector)
% - inclinations_indices (int32 vector)
% - azimuthal_indices (int32 vector)
%
% Usage:
% out = THDF_read_angular_grid_from_hdf5('path/to/file.h5')

% Auto-detect naming convention by checking which groups/datasets exist
infoRoot = h5info(h5file);
rootNames = {infoRoot.Groups.Name};
if isfield(infoRoot, 'Datasets') && ~isempty(infoRoot.Datasets)
    rootDsNames = {infoRoot.Datasets.Name};
else
    rootDsNames = {};
end

hasOutputField = any(contains(rootNames, '/output_field')) || any(contains(rootDsNames, 'output_field'));
hasEmergentField = any(contains(rootNames, '/emergent_field')) || any(contains(rootDsNames, 'emergent'));

if hasOutputField
    original_naming = true;
elseif hasEmergentField
    original_naming = false;
else
    error('Cannot detect HDF5 naming convention: neither /output_field nor /emergent_field found');
end

if hasOutputField
    original_naming = true;
elseif hasEmergentField
    original_naming = false;
else
    error('Cannot detect HDF5 naming convention: neither /output_field nor /emergent_field found');
end

if original_naming
    HDF5_NAMES = struct( ...
        'angular_grid', '/emergent_angular_grid', ...
        'output_frequencies', '/output_frequencies_grid', ...
        'N_directions', 'N_directions', ...
        'N_inclination_angles', 'N_inclination_angles', ...
        'N_azimuthal_angles', 'N_azimuthal_angles', ...
        'inclination_angles', 'inclination_angles', ...
        'azimuthal_angles', 'azimuthal_angles', ...
        'inclinations_indices', 'inclinations_indices', ...
        'azimuthal_indices', 'azimuthal_indices', ...
        'N_frequencies', 'N_frequencies', ...
        'frequencies', 'frequencies' ...
        );
else
    HDF5_NAMES = struct( ...
        'angular_grid', '/emergent_angular_grid', ...
        'output_frequencies', '/frequencies_grid', ...
        'N_directions', 'N_directions', ...
        'N_inclination_angles', 'N_inclination_angles', ...
        'N_azimuthal_angles', 'N_azimuthal_angles', ...
        'inclination_angles', 'inclination_angles', ...
        'azimuthal_angles', 'azimuthal_angles', ...
        'inclinations_indices', 'inclinations_indices', ...
        'azimuthal_indices', 'azimuthal_indices', ...
        'N_frequencies', 'N_frequencies', ...
        'frequencies', 'frequencies' ...
        );
end


if ~(ischar(h5file) || isstring(h5file))
    error('h5file must be a filename (string).');
end
h5file = char(h5file);
grp = HDF5_NAMES.angular_grid;

% Ensure group exists

try
    info = h5info(h5file, grp);
catch
    error(['Group ''' grp ''' not found']);
end

ds_names = {info.Datasets.Name};


% Read N_directions if present, otherwise infer from inclination_angles
if ismember(HDF5_NAMES.N_directions, ds_names)
    N_val = h5read(h5file, [grp '/' HDF5_NAMES.N_directions]);
    N = double(N_val(1));
else
    if ismember(HDF5_NAMES.inclination_angles, ds_names)
        incl_angles_tmp = h5read(h5file, [grp '/' HDF5_NAMES.inclination_angles]);
        N = size(incl_angles_tmp, 1);
    else
        error(['Neither ''' HDF5_NAMES.N_directions ''' nor ''' HDF5_NAMES.inclination_angles ''' found in group']);
    end
end


% Read N_inclination_angles and N_azimuthal_angles (expect scalars)
if ismember(HDF5_NAMES.N_inclination_angles, ds_names)
    N_incl_val = h5read(h5file, [grp '/' HDF5_NAMES.N_inclination_angles]);
    N_inclination_angles = int32(N_incl_val(1));
else
    error(['Dataset ''' HDF5_NAMES.N_inclination_angles ''' not found']);
end

if ismember(HDF5_NAMES.N_azimuthal_angles, ds_names)
    N_azim_val = h5read(h5file, [grp '/' HDF5_NAMES.N_azimuthal_angles]);
    N_azimuthal_angles = int32(N_azim_val(1));
else
    error(['Dataset ''' HDF5_NAMES.N_azimuthal_angles ''' not found']);
end


% Read arrays
inclination_angles = double(h5read(h5file, [grp '/' HDF5_NAMES.inclination_angles]));
azimuthal_angles = double(h5read(h5file, [grp '/' HDF5_NAMES.azimuthal_angles]));
inclinations_indices = int32(h5read(h5file, [grp '/' HDF5_NAMES.inclinations_indices]));
azimuthal_indices = int32(h5read(h5file, [grp '/' HDF5_NAMES.azimuthal_indices]));


% Try to read /output_frequencies_grid (optional)
freq_grp = HDF5_NAMES.output_frequencies;
try
    finfo = h5info(h5file, freq_grp);
    f_ds = {finfo.Datasets.Name};
    if ismember(HDF5_NAMES.N_frequencies, f_ds)
        Nf_val = h5read(h5file, [freq_grp '/' HDF5_NAMES.N_frequencies]);
        N_frequencies = int32(Nf_val(1));
    elseif ismember(HDF5_NAMES.frequencies, f_ds)
        tmpf = h5read(h5file, [freq_grp '/' HDF5_NAMES.frequencies]);
        N_frequencies = int32(size(tmpf, 1));
    else
        N_frequencies = int32(0);
    end
    if ismember(HDF5_NAMES.frequencies, f_ds)
        frequencies = double(h5read(h5file, [freq_grp '/' HDF5_NAMES.frequencies]));
    else
        frequencies = [];
    end
catch
    % frequencies group not present
    N_frequencies = int32(0);
    frequencies = [];
end

% Return struct
out = struct( ...
    'N_directions', int32(N), ...
    'inclination_angles', inclination_angles, ...
    'azimuthal_angles', azimuthal_angles, ...
    'inclinations_indices', inclinations_indices, ...
    'azimuthal_indices', azimuthal_indices, ...
    'N_inclination_angles', N_inclination_angles, ...
    'N_azimuthal_angles', N_azimuthal_angles, ...
    'N_frequencies', N_frequencies, ...
    'frequencies', frequencies ...
    );
end
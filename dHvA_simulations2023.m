% Initialize variables
band_numbers = [];       % To store the band numbers

% Check if the "program_runtime_settings" folder exists, and create it if necessary
settings_folder = 'program_runtime_settings';
if ~exist(settings_folder, 'dir')
    mkdir(settings_folder);
end

% Check if the user wants to use the same files as last time
use_last_files = input('Same files as last time? (y/n): ', 's');

% Load the file path list and band number list from the text file
settings_file = fullfile(settings_folder, 'last_run_files.txt');

if and(lower(use_last_files) == 'y', exist(settings_file, 'file'))
    % import each line in a cell
    saved_data = regexp(fileread(settings_file), '\r\n|\r|\n', 'split');
    material_name = saved_data{end - 1};
    band_numbers = str2num(saved_data{end});
    
    % initialize and populate the selected_files list
    selected_files = strings(1, numel(band_numbers));
    for i = 1:numel(band_numbers)
        selected_files(i) = saved_data{i};
    end
    % TODO Check if the file contents make sense and have been loaded correctly
    fprintf('Using the files from the last run.\n');
else
    % If the answer was yes, explain the problem
    if lower(use_last_files) == 'y'
        fprintf('The last_run_files file does not exist, searching for all SKEAF files...');
    end
    % Search for files in the directory tree
    start_directory = pwd;  % Get the current directory
    
    % recursively find files (gives a list of structs)
    files = dir(fullfile(start_directory, '**', 'results_freqvsangle.out'));
    % initialize the filepath list
    file_list = strings(1,length(files));

    for i = 1:length(files)
        file_list(i) = fullfile(files(i).folder, files(i).name);
    end
    
    % Display the found files with indices
    fprintf('Found %d SKEAF files:\n', numel(file_list));
    for i = 1:numel(file_list)
        fprintf(' - File %d: %s\n', i, file_list(i));
    end
    
    % Prompt the user to select files
    user_input = input('Which of the listed files should be read and combined? (give a list of integers, separated by space): ', 's');
    selected_indices = str2num(user_input);
    
    % Validate user input
    if isempty(selected_indices) || any(selected_indices < 1) || any(selected_indices > numel(file_list))
        error('Invalid input. Please enter valid integers within the range.');
        return;
    end
    
    % Create a list of selected file paths
    selected_files = file_list(selected_indices);
    fprintf('Selected files:\n');
    for i = 1:numel(selected_files)
        fprintf(' - File %d: %s\n', i, selected_files(i));
    end
    
    bands_match = false;
    % Ask for band numbers
    while ~bands_match
        band_input = input('What are the band numbers for these files? (give a list of integers, separated by space): ', 's');
        band_numbers = str2num(band_input);
        % check if the number of bands is the same as the number of files
        % selected. If not, make the user reintroduce the bands

        if numel(band_numbers) == numel(selected_indices)
            bands_match = true;
        else
            fprintf("The number of band numbers must match the number of selected files (%d).\n", numel(selected_indices));
        end
    end
    
    % ask for material name
    material_name = input("And finally, what is the name of the material? For plot titles. ", 's');

    % Save file paths and band numbers for future use
    file_data = sprintf('%s\n%s\n%s', strjoin(selected_files, '\n'), material_name, num2str(band_numbers));
    if ~exist('program_runtime_settings', 'dir')
        mkdir('program_runtime_settings');
    end

    fid = fopen('program_runtime_settings/last_run_files.txt', 'w');
    fprintf(fid, '%s', file_data);
    fclose(fid);
end

df = OrbitData();
df.set_material_name(material_name);
df.load_skeaf_files(selected_files, band_numbers);
% df.derivatives_calc(["freq_cos", "mass_cos"], 0.95);

% define the magnetic field range
Bmin = 5;
Bmax = 16;
BinvStep = 5e-5;
BinvRange = 1/Bmax: BinvStep: 1/Bmin;
BinvRange = BinvRange(end:-1:1);
Brange = 1./(BinvRange);

df.plot_freq_vs_angle_3d("mass_cos", true, df.material);
df.calculate_torque(Bmin, Bmax, BinvStep);
df.plot_torque_vs_field(Bmin, Bmax, BinvStep, 40, true);
df.calculate_torque_FFT(Bmin, Bmax, BinvStep, 1, 2);
df.plot_FFT(true, df.material, "freqMin", 1500, "freqMax", 3000, "angleMin", -0.1, "angleMax", 80);
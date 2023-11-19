classdef OrbitData < handle

    properties(SetAccess=private, GetAccess=public)
        material string
        data table %= table('Size', [0,9], 'VariableNames', ...
            %{'theta', 'phi', 'freq', 'freq_cos', 'mass', 'mass_cos', 'curv', 'dFdPhi', 'band'}, ...
            %'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'uint16'})
        torqueFFT
    end
    
    methods(Access=private)

    end

    methods(Access=public)
        function df = OrbitData()
            
        end

        function df = load_skeaf_files(df, filepath_list, band_numbers)
            
            % Initialize an empty cell array to store the data from each file
            data_cell = cell(length(filepath_list), 1);
        
            % Loop through each file in the list
            for i = 1:length(filepath_list)
                % Read data from the file, skipping the header
                skeaf_data = readmatrix(filepath_list{i}, 'FileType', 'text');
                
                % make mass negative if the orbit is hole-like
                skeaf_data(:, 4) = skeaf_data(:, 4) .* skeaf_data(:, 6);

                % get rid of orbit type and number of orbit columns
                skeaf_data = skeaf_data(:, 1:5);
        
                % Determine the number of rows in the data
                num_rows = size(skeaf_data, 1);
        
                % Create a vector with the corresponding band number for this file
                band_no = band_numbers(i) * ones(num_rows, 1);
        
                % Append the band number as a new column to the data
                skeaf_data = [band_no, skeaf_data];
        
                % Store the data in the cell array
                data_cell{i} = skeaf_data;
            end
        
            % Concatenate the data from all files into a single matrix
            data_cell = vertcat(data_cell{:});

            % order with increasing phi
            data_cell = sortrows(data_cell, 3);
            
            % assign data to the table
            df.data.band = data_cell(:, 1);
            df.data.theta = data_cell(:, 2);
            df.data.phi = data_cell(:, 3);
            df.data.freq = data_cell(:, 4) .* 1000;
            df.data.freq_cos = df.data.freq .* cos(pi .* df.data.phi ./ 180);
            df.data.mass = data_cell(:, 5);
            df.data.mass_cos = df.data.mass .* cos(pi .* df.data.phi ./ 180);
            df.data.curv = data_cell(:, 6);
        end
        
        function plot_freq_vs_angle(df, cosine, freqCosMin, freqCosMax, angleMin, angleMax, varargin)
            % check if inputs are fine
            if ~isnumeric(freqCosMin) || ~isnumeric(freqCosMax) || ~isnumeric(angleMin) || ~isnumeric(angleMax)
                error("plot_freq_vs_angle did not receive numeric inputs");
            end

            if ~islogical(cosine)
                error("plot_freq_vs_angle cosine argument received is not a logical value");
            end
            
            % prepare mask
            mask = (df.data.phi > angleMin) & (df.data.phi < angleMax) &...
            (df.data.freq_cos > freqCosMin) & (df.data.freq_cos < freqCosMax);

            % if Fcos plot, make sure the ydata is correct and apply mask
            if cosine
                ydata = df.data.freq_cos(mask);
                freq_label = "Fcos\phi (T)";
            else
                ydata = df.data.freq(mask);
                freq_label = "Frequency (T)";
            end
            
            xdata = df.data.phi(mask);
            
            % plot away
            figure;
            gscatter(xdata, ydata, df.data.band(mask), [], "+");
            xlabel("Angle \phi (degrees)");
            ylabel(freq_label);

            % lazy way to set title
            title(varargin{1});
        end
        
        function plot_freqcos_vs_angle_mass(df, freqCosMin, freqCosMax, angleMin, angleMax, massMin, massMax, varargin)
            % check if inputs are fine
            if ~isnumeric(freqCosMin) || ~isnumeric(freqCosMax) ||...
                    ~isnumeric(angleMin) || ~isnumeric(angleMax) ||...
                    ~isnumeric(massMin) || ~isnumeric(massMax)
                error("plot_freq_vs_angle did not receive numeric inputs");
            end

            % prepare mask
            mask = (df.data.phi > angleMin) & (df.data.phi < angleMax) &...
            (df.data.freq_cos > freqCosMin) & (df.data.freq_cos < freqCosMax) &...
            (df.data.mass > massMin) & (df.data.mass < massMax);
            
            % apply mask
            xdata = df.data.phi(mask);
            ydata = df.data.freq_cos(mask);
            zdata = df.data.mass(mask);
            
            % plot in 3D, but set view as 2D
            figure;
            scatter3(xdata, ydata, zdata, [], zdata);
            view(2);
            colorbar;
            xlabel("Angle \phi (degrees)");
            ylabel("Fcos\phi (T)");
            zlabel("Mass (electron mass units)");
            third_axis_colorbar = colorbar(gca);
            third_axis_colorbar.Label.String = "Mass (electron mass units)";
            title(varargin{1});
        end

        function plot_torque_vs_field(df, torque_of_orbits, Brange, angleMin, angleMax, torque_vs_angle_title)
            if ~isnumeric(torque_of_orbits) || ~isnumeric(Brange) || ~isnumeric(angleMin)...
                    || ~isnumeric(angleMax) || ~isstring(torque_vs_angle_title)
                error("plot_torque_vs_field parameters are not the right type")
            end
            
            angle_mask = (df.data.phi > angleMin - 0.00001) & (df.data.phi < angleMax + 0.00001);
            
            torque_masked = torque_of_orbits(angle_mask, :);
            torque_response_at_angle = sum(torque_masked, 1);
            
            plot(Brange, torque_response_at_angle);
            title(torque_vs_angle_title);
            xlabel("Magnetic field (T)");
            ylabel("Torque (arbitrary units)");
            
        end
        
        function df = derivatives_calc(df, features_used, k_dist_perc)
            % Written by Alexandru Dobra, September 2023
            %
            % Inputs:
            % df: a table containing data from the SKEAF simulation
            % (df has to have columns named phi and freq, which is the polar angle in my case, along
            % which the derivative is calculated)
            % features_used: a cell array containing the names of the columns one wishes to group by
            % (in my case, mass and freq or mass_cos and freq_cos)
            % k_dist_perc: the percentile distance to the 3rd nearest neighbor -
            % this is used to determine the cutoff distance, I use 0.95
            %
            % Output: the same df, with an extra column called dFdPhi in which the derivatives
            % are stored
        
            % Extract the list of unique angles
            angleList = unique(df.data.phi);
            angleNumber = numel(angleList);
        
            % Initialize the derivatives column
            df.data.dFdPhi = cell(height(df.data), 1);
        
            for i = 2:(angleNumber - 1)
                % Get the datapoints of the current slice
                angleMask = (df.data.phi > (angleList(i - 1) - 0.01)) & (df.data.phi < (angleList(i + 1) + 0.01));
                % Group them based on proximity
                slice = df.data{angleMask, features_used};
                start_idx = find(angleMask, 1);
                groups = groupEmUp(slice, k_dist_perc);
        
                valid_groups = {};
                for j = 1:length(groups)
                    group = groups{j} + start_idx - 1;
                    % Remove groups that contain all datapoints at the same angle
                    if abs(max(df.data.phi(group)) - min(df.data.phi(group))) < 0.00001
                        continue;
                    end
                    valid_groups = [valid_groups, group];
                end
        
                groups = valid_groups;
        
                % Initialize derivatives list
                derivatives = cell(1, length(groups));
        
                for j = 1:length(groups)
                    group = groups{j};
        
                    % Approximate the derivative by fitting to a line
                    xdata = df.data.phi(group);
                    ydata = df.data.freq(group);
                    [a, ~] = polyfit(xdata, ydata, 1);
                    derivatives{j} = a(1);
                end
        
                % Append derivative values to each datapoint of the group in the big table
                for j = 1:length(groups)
                    group = groups{j};
                    df.data.dFdPhi(group) = cellfun(@(x) [x, derivatives{j}], df.data.dFdPhi(group), 'UniformOutput', false);
                end
            end
        
            % Remove datapoints that have no dFdPhi values
            for i = 1:height(df.data)
                if isempty(df.data.dFdPhi{i})
                    df.data.dFdPhi{i} = NaN;
                end
            end
        
            % Remove rows with NaN in dFdPhi column
            df.data = df.data(~cellfun(@(x) all(isnan(x)), df.data.dFdPhi), :);
        
            % Average the derivative values found
            for i = 1:height(df.data)
                df.data.dFdPhi{i} = mean(df.data.dFdPhi{i});
            end
        
            df.data.dFdPhi = cell2mat(df.data.dFdPhi);
        end
        
        function torque = calculate_torque(df, Brange, pmax)
            if ~isnumeric(Brange) || ~isnumeric(pmax)
                error("calculate_torque parameters are not the right type")
            end

            % load relevant constants from a file
            load physicalConstants-SI.mat e hbar me kb
            % define temperature at 2K
            temperature = 2;
            % set the mass reduction factor to 0.1
            mm = 0.1 .* me;
            % set the spin-splitting factor
            g = 2;
            % set the mean free path (in Angstrom)
            mean_free_path = 300;
            % set berry phase
            berry_phase = 0;

            slice = df.data(:, :);
            
            % initialize torque with zeros in the correct dimensions
            torque = zeros(height(slice), length(Brange));

            for p = 1:pmax
                chi = 2.*pi.^2.*p.*kb/(e.*hbar).*temperature.*mm...
                    .* slice.mass * (Brange.^(-1));
                % define the damping factors
                RT = chi ./ sinh(chi);
                RS = cos(.5.*pi.*p.*g.* slice.mass) * ones(1, length(Brange));
                RD = exp(-1140.*p./mean_free_path .* (slice.freq.^0.5) * (Brange.^(-1)));
                % add the current harmonic's influence on the torque
                torque = torque + p.^(-1.5) .* RT.*RS.*RD .*...
                    sin(2.*pi.*p.* ...
                        (slice.freq * (Brange.^(-1)) - 0.5 + berry_phase/2/pi +...
                            1/8.*sign(slice.curv) * ones(1, length(Brange))...
                        )...
                    );
            end
            
            % multiply with the factor independent of p
            torque = ((abs(slice.curv).^(-.5)) * (Brange.^(3/2))).*...
                (slice.dFdPhi * ones(1, length(Brange))) .* ...
                (slice.mass.^(-1) * ones(1, length(Brange))) .* torque;
        end

        function calculate_torque_FFT(df, Bmin, Bmax, BinvStep, pmax)
            % make the B-field uniform in 1/B as that gives uniform
            % frequency when taking the FFT
            BinvRange = 1/Bmax: BinvStep: 1/Bmin;
            BinvRange = BinvRange(end:-1:1);
            Brange = 1./(BinvRange);
            
            % calculate torque for all angles
            angleMin = 0;
            angleMax = 90;
            torque = calculate_torque(df, Brange, angleMin, angleMax, pmax);
            df.torqueFFT = abs(fft(torque, [], 2))./abs(BinvRange(end) - BinvRange(1));
            
        end


    end
end

% Torque response for CsV_3Sb_5 no spin-orbit bands 67 & 68, at \theta=30 deg 
% Magnetic field (T)
% Torque (arbitrary units)

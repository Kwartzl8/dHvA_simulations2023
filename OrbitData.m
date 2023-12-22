classdef OrbitData < handle

    properties(SetAccess=private, GetAccess=public)
        material string
        data table
        torque
        FFT_freq_range
        torqueFFT
    end
    
    methods(Access=private)
        function w = hannwin(df, L)
        n = 1:1:L;
        w=(0.5*(1-cos(2*pi*n/L)));
        end
        
        function w = gausswin(df, L)
        N = L-1;
        n = -N/2:1:N/2;
        alpha = 1.5;
        win1 = exp(-0.5*(2*alpha.*(n)/N).^2);
        w =win1;
        end
        
        function w = hammingwin(df, L)
        n = 1:1:L;
        w=(25/46*(1-cos(2*pi*n/L)));
        end
        
        function w = blackmanwin(df, L)
        % using exact blackman
        n = 1:1:L;
        a0 = 7938/18608;
        a1 = 9240/18608;
        a2 = 1430/18608;
        w = a0 - a1*cos(2*pi*n/L) + a2*cos(4*pi*n/L); 
        end
        
        function w = welchwin(df, L)
        n = 1:1:L;
        w=1 - ((n-L/2)/(L/2)).^2;
        end
    end

    methods(Access=public)
        function df = OrbitData()
            
        end

        function set_material_name(df, material_name)
            df.material = material_name;
        end

        function load_skeaf_files(df, filepath_list, band_numbers)
            
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

        function plot_freq_vs_angle_3d(df, zdata_name, freq_cos_flag, plot_title, limits)
            arguments
                df (1, 1) OrbitData
                zdata_name (1,:) char {mustBeMember(zdata_name, ...
                    {'theta', 'phi', 'freq', 'freq_cos', 'mass', ...
                    'mass_cos', 'curv', 'dFdPhi', 'band'})} = 'mass'
                freq_cos_flag (1, 1) logical = true
                plot_title (1, 1) string {} = zdata_name + " plot"
                limits.freqMin (1, 1) double = 0
                limits.freqMax (1, 1) double = 1e5
                limits.angleMin (1, 1) double = -0.1
                limits.angleMax (1, 1) double = 90.1
                limits.zMin (1, 1) double = -Inf
                limits.zMax (1, 1) double = Inf
            end

            mask = (df.data.phi > limits.angleMin) & (df.data.phi < limits.angleMax) &...
                (df.data.freq_cos > limits.freqMin) & (df.data.freq_cos < limits.freqMax) &...
                table2array(df.data(:, zdata_name) > limits.zMin) & table2array(df.data(:, zdata_name) < limits.zMax);

            xdata = df.data.phi(mask);

            if freq_cos_flag
                ydata = df.data.freq_cos(mask);
            else
                ydata = df.data.freq(mask);
            end

            zdata = table2array(df.data(mask, zdata_name));
            
            figure;
            scatter3(xdata, ydata, zdata, [], zdata);
            xlabel("Angle $\phi$ (degrees)", 'Interpreter','latex');

            if freq_cos_flag
                ylabel("$F\cos{\phi}$ (T)",'Interpreter','latex');
            else
                ylabel("Frequency $F$ (T)",'Interpreter','latex');
            end

            zlabel(zdata_name);
            cb = colorbar();
            ylabel(cb, zdata_name);
            title(plot_title, "Interpreter", "latex");
            set(gca,'TickLabelInterpreter','latex');
        end

        function plot_torque_vs_field(df, Bmin, Bmax, BinvStep, angle, invert_B, plot_title)
            arguments
                df OrbitData
                Bmin double
                Bmax double
                BinvStep double
                angle double
                invert_B logical = false
                plot_title string = df.material + " torque vs field"
            end
            
            if isempty(df.torque)
                error("Torque data is missing. Call calculate_torque on this object.");
            end
            
            % prepare x-axis
            BinvRange = 1/Bmax: BinvStep: 1/Bmin;
            Brange = BinvRange(end:-1:1);
            if ~invert_B
                Brange = 1./(Brange);
            end
            
            % prepare y-axis
            distinct_angles = unique(df.data.phi);
            [~, closest_angle_index] = min(abs(distinct_angles - angle));
            closest_angle = distinct_angles(closest_angle_index);
            
            torque_at_given_angle = df.torque(closest_angle_index, :);
            
            % plot
            figure;
            plot(Brange, torque_at_given_angle, "DisplayName", "angle " + num2str(closest_angle, 3));
            title(plot_title, 'Interpreter', 'latex');
            legend();

            if invert_B
                xlabel_string = "Inverse magnetic field $1/B$ (1/T)";
            else
                xlabel_string = "Magnetic field $B$ (T)";
            end
            xlabel(xlabel_string, 'Interpreter', 'latex');
            ylabel("Torque $\tau$ (arbitrary units)", 'Interpreter', 'latex');
            set(gca,'TickLabelInterpreter','latex');
        end

        function plot_FFT(df, freq_cos_flag, plot_title, scaling, save_fig, filename, limits)
            arguments
                df (1, 1) OrbitData
                freq_cos_flag (1, 1) logical = true
                plot_title (1, 1) string {} = zdata_name + " plot"
                scaling (1, :) char {mustBeMember(scaling, {'none', 'sqrt', 'log'})} = 'sqrt'
                save_fig (1, 1) logical = false
                filename (1, :) char = 'torque_FFT_plot'
                limits.freqMin (1, 1) double = 0
                limits.freqMax (1, 1) double = 1e5
                limits.angleMin (1, 1) double = -0.1
                limits.angleMax (1, 1) double = 90.1
            end

            % has the FFT been calculated?
            if isempty(df.torqueFFT)
                error("The torque has not been calculated for this OrbitData object. Call calculate_torque_FFT on the object first.");
            end

            distinct_angles = unique(df.data.phi);

            [freq_mesh, angle_mesh] = meshgrid(df.FFT_freq_range, distinct_angles);
            [number_of_distinct_angles, number_of_distinct_frequencies] = size(freq_mesh);

            if freq_cos_flag
                freq_mesh = freq_mesh .* (cos(pi .* distinct_angles ./ 180) * ones(1, number_of_distinct_frequencies));
            end
            
            zdata = df.torqueFFT;
            switch scaling
                case "none"
                    "nothing";
                case "sqrt"
                    zdata = sqrt(zdata);
                case "log"
                    zdata = log10(zdata + 1);
            end

            figure;
            waterfall(angle_mesh, freq_mesh, zdata);
            xlim([limits.angleMin, limits.angleMax]);
            ylim([limits.freqMin, limits.freqMax]);

            xlabel("Angle $\phi$ (degrees)", 'Interpreter','latex');

            if freq_cos_flag
                ylabel("$F\cos{\phi}$ (T)",'Interpreter','latex');
            else
                ylabel("Frequency $F$ (T)",'Interpreter','latex');
            end
            
            zlabel_string = "Torque $\tau$ (arbitraty units)";
            zlabel(zlabel_string, 'Interpreter', 'latex');
            cb = colorbar();
            ylabel(cb, zlabel_string, 'Interpreter', 'latex');
            title(plot_title, 'Interpreter', 'latex');
            set(gca,'TickLabelInterpreter','latex');

            if save_fig
                % if no filename given, add the time to it for a lower
                % likelihood of overwriting
                if filename == "torque_FFT_plot"
                    filename = filename + datetime('now', 'Format', 'HH_mm_ss');
                end
                
                if ~exist('saved_figures', 'dir')
                    mkdir('saved_figures');
                end
                savefig("./saved_figures/" + filename);
            end
        end

        function derivatives_calc(df, features_used, k_dist_perc)
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
        
        function calculate_torque(df, Bmin, Bmax, BinvStep, pmax, parameters)
        arguments
            df OrbitData
            Bmin (1,1) double
            Bmax (1,1) double
            BinvStep (1,1) double
            pmax (1,1) double = 1
            parameters.Temperature double = 0.45
            parameters.MassReductionFactor double = 0.1
            parameters.MeanFreePath double = 300
            parameters.BerryPhase double = 0
        end
            
            % check if the derivative has been calculated
            if ~strcmp(df.data.Properties.VariableNames, "dFdPhi")
                warning("Derivative values are missing. Calculating with default parameters...")
                df.derivatives_calc(["freq_cos", "mass_cos"], 0.95);
            end

            % load relevant constants from a file
            load physicalConstants-SI.mat e hbar me kb
            
            temperature = parameters.Temperature;
            % set the mass reduction factor
            mm = parameters.MassReductionFactor .* me;
            % set the spin-splitting factor
            g = 2;
            % set the mean free path (in Angstrom)
            mean_free_path = parameters.MeanFreePath;
            % set berry phase
            berry_phase = parameters.BerryPhase;
            
            BinvRange = 1/Bmax: BinvStep: 1/Bmin;
            BinvRange = BinvRange(end:-1:1);
            Brange = 1./(BinvRange);
            slice = df.data(:, :);
            
            % initialize torque with zeros in the correct dimensions
            torque_of_orbits = zeros(height(slice), length(Brange));

            for p = 1:pmax
                chi = 2.*pi.^2.*p.*kb/(e.*hbar).*temperature.*mm...
                    .* slice.mass * (Brange.^(-1));
                % define the damping factors
                RT = chi ./ sinh(chi);
                RS = cos(.5.*pi.*p.*g.* slice.mass .* mm/me) * ones(1, length(Brange));
                % RS = ones(height(slice), length(Brange));
                RD = exp(-1140.*p.^0.5./mean_free_path .* (slice.freq.^0.5) * (Brange.^(-1)));
                % add the current harmonic's influence on the torque
                torque_of_orbits = torque_of_orbits + p.^(-1.5) .* RT.*RS.*RD .*...
                    sin(2.*pi.*p.* ...
                        (slice.freq * (Brange.^(-1)) - 0.5 + berry_phase/2/pi +...
                            1/8/p.*sign(slice.curv) * ones(1, length(Brange))...
                        )...
                    );
            end
            
            % multiply with the factor independent of p
            torque_of_orbits = ((abs(slice.curv).^(-.5)) * (Brange.^(3/2))).*...
                (slice.dFdPhi * ones(1, length(Brange))) .* ...
                (slice.mass.^(-1) * ones(1, length(Brange))) .* torque_of_orbits;

            % group the values based on angle
            sum_along_dimension1 = @(matrix) sum(matrix, 1);
            angle_groups = findgroups(df.data.phi);
            df.torque = splitapply(sum_along_dimension1, torque_of_orbits, angle_groups);
        end

        function calculate_torque_FFT(df, Bmin, Bmax, BinvStep, pmax, win)
            % make the B-field uniform in 1/B as that gives uniform
            % frequency when taking the FFT
            BinvRange = 1/Bmax: BinvStep: 1/Bmin;
            BinvRange = BinvRange(end:-1:1);
            Brange = 1./(BinvRange);
            Brange_length = length(Brange);
            
            df.FFT_freq_range = 1/(1/Brange(1) - 1/Brange(end))*(0:round(Brange_length/2));

            % torque has dimensions of len(unique(df.data.phi)) x len(BinvRange)
            if isempty(df.torque)
                warning("Torque data is missing. Calculating torque with the default parameters...");
                df.calculate_torque(Bmin, Bmax, BinvStep, pmax);
            end
            
            unique_angles = unique(df.data.phi);
            % windowing
            %   1=None 2=Hanning 3=Gauss 4=Hamming 5=Blackman 6=welch
            if win == 2
                w = df.hannwin(Brange_length);
            elseif win == 3
                w = df.gausswin(Brange_length);
            elseif win == 4
                w = df.hammingwin(Brange_length);
            elseif win == 5
                w = df.blackmanwin(Brange_length); 
            elseif win == 6
                w = df.welchwin(Brange_length);
            else
                w = ones(1, Brange_length);     % no windowing
            end
            % scale to matrix
            [~, w_mesh] = meshgrid(unique_angles, w);

            torque_windowed = w_mesh.' .* df.torque;

            df.torqueFFT = abs(fft(torque_windowed, [], 2))./abs(BinvRange(end) - BinvRange(1));
            df.torqueFFT = df.torqueFFT(:, 1:round(Brange_length/2)+1);
        end

    end
end

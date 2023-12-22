# OrbitData Class Documentation
## Class Definition

The `OrbitData` class in MATLAB is designed to handle and analyze data related to dHvA orbits. It provides functionality for loading data, calculating torque, performing Fourier transforms, and visualizing results.

### Properties

- `material` (string): The material name associated with the data.
- `data` (table): A table containing orbit data, including columns for band, theta, phi, frequency, mass, curvature, and additional derived columns like frequency times the cosine of phi.
- `torque`: Placeholder for torque data calculated using the `calculate_torque` method.
- `FFT_freq_range`: Frequency range array for the Fourier transform, filled in the `calculate_torque_FFT` method.
- `torqueFFT`: Placeholder for Fourier-transformed torque data calculated using the `calculate_torque_FFT` method.

### Private Methods

#### Windowing Functions
These are used when calculating the Fourier transform in `calculate_torque_FFT`.
- `hannwin(df, L)`: Hanning window function.
- `gausswin(df, L)`: Gaussian window function.
- `hammingwin(df, L)`: Hamming window function.
- `blackmanwin(df, L)`: Blackman window function.
- `welchwin(df, L)`: Welch window function.

### Public Methods

#### Constructor

- `df = OrbitData()`: Creates an instance of the `OrbitData` class.

#### Data Loading and Manipulation

- `set_material_name(df, material_name)`: Sets the material name.
- `load_skeaf_files(df, filepath_list, band_numbers)`: Loads SKEAF simulation files, processes data, and populates the `data` table.

#### Visualization

- `plot_freq_vs_angle_3d(df, zdata_name, freq_cos_flag, plot_title, limits)`: Creates a 3D scatter plot of frequency vs angle, the z-axis variable can be chosen by the user out of the columns of the `df.data` table.
- `plot_torque_vs_field(df, Bmin, Bmax, BinvStep, angle, invert_B, plot_title)`: Plots torque vs magnetic field at the angle closest to the one given. It can also plot torque vs inverse magnetic field (in which the periodicity of dHvA oscillations is apparent) if the option `invert_B` is set to `true`.
- `plot_FFT(df, freq_cos_flag, plot_title, scaling, save_fig, filename, limits)`: Plots the Fourier transform of torque.

#### Derivative Calculation

- `derivatives_calc(df, features_used, k_dist_perc)`: Calculates derivatives of frequency with respect to the polar angle phi. I recommend setting the `features_used` parameter to `["freq_cos", "mass_cos"]` and `k_dist_perc` to `0.95`.

#### Torque Calculation

- `calculate_torque(df, Bmin, Bmax, BinvStep, pmax, parameters)`: Calculates torque based on specified parameters. `pmax` is the maximum number of harmonics - default is `1`. Other parameters can be changed including temperature, the mass reduction factor and the mean free path. The result is saved in `df.torque` and has dimensions of `length(unique(df.data.phi)) x length(BinvRange)`, meaning that the torque was calculated for each angle and each magnetic field value. The torque is calculated using the following formula:
$$
    \tau=\left(\frac{e}{\hbar}\right)^{\frac{3}{2}}\left(\frac{e \hbar}{m^* c}\right) \frac{B^{\frac{3}{2}} V}{2^{\frac{1}{2}} \pi^{\frac{3}{2}}\left|\frac{d^2A}{dk^2_{\parallel}}\right|^{\frac{1}{2}}} \frac{\partial F}{\partial \theta} \sum_{p=1}^{\infty} p^{-\frac{3}{2}} R_S R_D R_T \sin \left(2 \pi p\left(\frac{F}{B}-\frac{1}{2} + \frac{1}{8}\text{sgn}{\frac{d^2A}{dk^2_{\parallel}}} +\frac{\phi_B}{2 \pi}\right)\right)
$$
	where the symbols $e$, $\hbar$ and $c$ have their usual meaning, $B$ is the magnetic field, $V$ is the sample volume, $\frac{d^2A}{dk^2_{\parallel}}$ is the curvature mentioned before, $F$ is the orbit frequency, $\phi_B$ is the Berry phase and the damping factors $R_S$, $R_D$ and $R_T$ are given by:
$$\begin{aligned}
    R_T & =\frac{\chi}{\sinh (\chi)} \quad \chi=\frac{2 \pi^2 p k_B T m^*}{e \hbar B} \\ R_S & =\cos \left(1 / 2 p \pi g m^* / m_e\right) \\ R_D & =\exp \left(\frac{-2 \pi^2 p k_B x}{e\hbar B}m^*c\right)=\exp \left(-1140 p \sqrt{F} / (l B)\right)\end{aligned}$$
	where $g$ is the spin-splitting factor, $x = \hbar / 2 \pi k_B \tau$ (with $\tau$ the electron scattering time) is the Dingle temperature (depends on the level of impurities in the material), and $l$ is the mean free path (in Angstrom). Keep in mind that the code omits some of the fundamental constants and the volume of the sample for clarity. The result is in arbitrary units of torque.
- `calculate_torque_FFT(df, Bmin, Bmax, BinvStep, pmax, win)`: Calculates the Fourier transform of torque in the B-field direction. The result is stored in `df.torqueFFT` and has dimensions of `length(unique(df.data.phi)) x length(df.FFT_freq_range)`
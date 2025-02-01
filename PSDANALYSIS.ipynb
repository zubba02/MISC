import numpy as np
import scipy.io
from scipy import spatial
import matplotlib.pyplot as plt

mat = scipy.io.loadmat('TEST_ne_e.mat')

x_vals = mat['Xp'].ravel()
y_vals = mat['Yp'].ravel()
coords  = [[x_vals[i], y_vals[i]] for i in range(0, np.shape(x_vals)[0])]
tree = spatial.KDTree(coords)

coordslist = np.genfromtxt('PROF_1.csv', delimiter=',')

prof = 'PROFILE_1'

Hs_all = []
Tp_all = []

for z in coordslist:
    print (z)
    extract_coords = z[1],z[2]
    val = tree.query([(extract_coords[0], extract_coords[1])])[1]

    all_keys = mat.keys()
    all_keys_watlev = [i for i in all_keys if i[0:3] == 'Wat']

    waterlevels = []

    for j in all_keys_watlev:
        vals = mat['{}'.format(j)].ravel()[val]
        waterlevels.append(vals)

    plt.figure(figsize=(15, 6))
    time_vals = np.arange(0, 2162.0, 2)
    plt.plot(time_vals,waterlevels)
    plt.title("Water Elevation at Location {} {}".format(z[0], prof))
    plt.xlabel("Time (s)")
    plt.ylabel("Water Elevation (m)")
    plt.savefig('Watlev_at Location {} {}.png'.format(z[0], prof), dpi=300)

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import welch
    from scipy.ndimage import gaussian_filter1d
    
    Times = np.array(time_vals)
    Watlev = np.array(waterlevels).flatten()  # Flatten to 1D array    

    sampling_interval = np.mean(np.diff(Times))  # Average time step
    sampling_frequency = 1 / sampling_interval  # Sampling frequency (Hz)
    
    frequencies, psd = welch(Watlev, fs=sampling_frequency, nperseg=1024)

    mask = (frequencies >= 0) & (frequencies <= 0.5)
    frequencies = frequencies[mask]
    psd = psd[mask]

    psd_normalized = psd / np.max(psd)

    psd_smooth = gaussian_filter1d(psd_normalized, sigma=2)  # Adjust sigma for smoothing strength

    plt.figure(figsize=(15, 6))
    plt.semilogy(frequencies, psd_smooth, label="Wave Spectrum (PSD)", color="blue")
    plt.title("Wave Spectrum of Watlev Time Series at Location {} PROF_1".format(z[0]))
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Spectral Density (m²/Hz)")
    plt.grid(True, which="both", ls="--")
    plt.xlim(0, 0.25)  # Set x-axis limits
    plt.ylim(0, 1)    # Set y-axis limits
    plt.legend()
    plt.savefig('Wavespectram_at Location_{}_{}.png'.format(z[0], prof), dpi=300)
    #plt.show()

    std_watlev = np.std(Watlev)
    Hs = 4 * std_watlev
    
    peak_frequency_index = np.argmax(psd_smooth)  # Index of maximum PSD value
    peak_frequency = frequencies[peak_frequency_index]  # Peak frequency
    Tp = 1 / peak_frequency  # Peak period

    Hs_all.append(Hs)
    Tp_all.append(Tp)

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import welch
    from scipy.ndimage import gaussian_filter1d

    # Step 1: Define the Times and Watlev arrays
    Times = np.array(time_vals)
    
    Watlev = np.array(waterlevels).flatten()  # Flatten to 1D array
    
    # Step 2: Compute sampling frequency
    sampling_interval = np.mean(np.diff(Times))  # Average time step
    sampling_frequency = 1 / sampling_interval  # Sampling frequency (Hz)
    
    # Step 3: Compute the Power Spectral Density (PSD) using Welch's method
    frequencies, psd = welch(Watlev, fs=sampling_frequency, nperseg=1024)
    
    # Step 4: Limit the frequency range to 0–0.5 Hz
    mask = (frequencies >= 0) & (frequencies <= 0.5)
    frequencies = frequencies[mask]
    psd = psd[mask]
    
    # Step 5: Normalize the PSD to the range 0–1
    psd_normalized = psd / np.max(psd)
    
    # Step 6: Smoothen the PSD using a Gaussian filter
    psd_smooth = gaussian_filter1d(psd_normalized, sigma=2)  # Adjust sigma for smoothing strength
    
    # Step 7: Define frequency ranges for wave types
    frequency_ranges = {
        "Capillary Waves": (0.5, np.inf),
        "Wind Waves": (0.1, 0.5),
        "Swell Waves": (0.05, 0.1),
        "Infragravity Waves": (0.002, 0.05),
        "Far Gravity Waves": (0.0, 0.002)
    }
    
    # Step 8: Compute energy in each frequency band
    total_energy = np.trapz(psd_smooth, frequencies)  # Total energy (area under the curve)
    energy_percentages = {}
    
    for wave_type, (f_min, f_max) in frequency_ranges.items():
        mask = (frequencies >= f_min) & (frequencies <= f_max)
        energy_in_band = np.trapz(psd_smooth[mask], frequencies[mask])  # Energy in the band
        percentage = (energy_in_band / total_energy) * 100  # Percentage of total energy
        energy_percentages[wave_type] = percentage

    WAVETYPES = []
    PERC = []
    # Step 9: Print results
    print("Energy Distribution Across Wave Types:")
    for wave_type, percentage in energy_percentages.items():
        #print(f"{wave_type}: {percentage:.2f}%")
        WAVETYPES.append(wave_type)
        PERC.append(percentage)
    print (WAVETYPES)
    
    # Step 10: Plot the wave spectrum with frequency bands
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, psd_smooth, label="Smoothed Wave Spectrum", color="blue")
    
    # Highlight frequency bands
    colors = ["red", "orange", "green", "purple", "brown"]
    for (wave_type, (f_min, f_max)), color in zip(frequency_ranges.items(), colors):
        if f_max > frequencies[-1]:  # Skip capillary waves if out of range
            continue
        plt.axvspan(f_min, f_max, color=color, alpha=0.2, label=f"{wave_type}")   

    plt.title('Wave Spectrum at Location_{}_{}.png'.format(z[0], prof))
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Normalized Spectral Density")
    plt.xlim(0, 0.5)  # Set x-axis limits
    plt.ylim(0, 1)    # Set y-axis limits
    plt.grid(True, which="both", ls="--")
    plt.text(0.35, 0.75, 'Capillary Waves: {0:.2f}%'.format(PERC[0]))
    plt.text(0.35, 0.7, 'Wind Waves: {0:.2f}%'.format(PERC[1]))
    plt.text(0.35, 0.65, 'Swell Waves: {0:.2f}%'.format(PERC[2]))
    plt.text(0.35, 0.60, 'Infragravity Waves: {0:.2f}%'.format(PERC[3]))
    plt.text(0.35, 0.55, 'Far Gravity Waves: {0:.2f}%'.format(PERC[4]))
    plt.text(0.35, 0.5, 'Hs: {0:.2f}m'.format(Hs))
    plt.text(0.35, 0.45, 'Tp: {0:.2f}s'.format(Tp))
    plt.legend()
    plt.savefig('Wavespectrum_FB_Location_{}_{}.png'.format(z[0], prof), dpi=300)
    #plt.show()


np.savetxt('Hs.csv', Hs_all, delimiter=',')
np.savetxt('Tp.csv', Tp_all, delimiter=',')

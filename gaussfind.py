import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define multimodal Gaussian function
def multimodal_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        a = params[i]
        mu = params[i+1]
        sigma = params[i+2]
        y += a * np.exp(-((x - mu)**2) / (2 * sigma**2))
    return y

# Fit function
def fit_multimodal_gaussian(x_data, y_data, n_modes, plot=True):
    # Initial guesses: evenly spaced means, all amplitudes ~max(y), all sigmas ~10
    a_guess = max(y_data)
    mu_guess = np.linspace(min(x_data), max(x_data), n_modes)
    sigma_guess = np.full(n_modes, 10)
    
    initial_params = []
    for i in range(n_modes):
        initial_params += [a_guess, mu_guess[i], sigma_guess[i]]

    # Fit
    popt, _ = curve_fit(multimodal_gaussian, x_data, y_data, p0=initial_params)

    # Optional plotting
    if plot:
        x_fit = np.linspace(min(x_data), max(x_data), 1000)
        y_fit = multimodal_gaussian(x_fit, *popt)

        plt.plot(x_data, y_data, 'bo', label='Data')
        plt.plot(x_fit, y_fit, 'r-', label='Fit')

        # Plot each individual Gaussian component
        for i in range(n_modes):
            a = popt[i*3]
            mu = popt[i*3+1]
            sigma = popt[i*3+2]
            y_component = a * np.exp(-((x_fit - mu)**2) / (2 * sigma**2))
            plt.plot(x_fit, y_component, '--', label=f'Component {i+1}')
        
        plt.legend()
        plt.xlabel('λ')
        plt.ylabel('Absorption')
        plt.title(f'Gaussian Fit with {n_modes} Modes')
        plt.show()

    # Return structured results
    components = []
    for i in range(n_modes):
        components.append({
            'amplitude': popt[i*3],
            'mean': popt[i*3+1],
            'width': popt[i*3+2]
        })
    return components
lambda_vals = np.array([
455.328,
443.443,
448.77,
435.656,
430.738,
427.459,
422.131,
414.344,
409.426,
404.918,
399.18,
391.393,
382.787,
373.77,
366.803,
360.656,
355.738,
350.82,
462.705,
465.574,
468.443,
469.672,
470.492,
472.541,
475,
477.869,
479.098,
481.557,
488.934,
484.426,
497.951,
343.852,




])
absorption = np.array([
0.484158,
0.444554,
0.472277,
0.387129,
0.332673,
0.30297,
0.250495,
0.193069,
0.168317,
0.153465,
0.155446,
0.176238,
0.194059,
0.177228,
0.154455,
0.128713,
0.10198,
0.080198,
0.455446,
0.425743,
0.369307,
0.339604,
0.309901,
0.253465,
0.19604,
0.137624,
0.108911,
0.080198,
0.0247525,
0.049505,
0.00792079,
0.0544554,




])
components = fit_multimodal_gaussian(x_data=lambda_vals, y_data=absorption, n_modes=2)
for i, comp in enumerate(components):
    print(f"Mode {i+1}: height={comp['amplitude']:.3f}, mean={comp['mean']:.2f}, width={comp['width']:.2f}")
# Example usage: use plot digitizer data to get absorption parameters.
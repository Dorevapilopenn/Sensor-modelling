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
453.02,
433.52,
380.787,
401.677,
422.337,
360.442,
369.839,
392.179,
464.529,
472.003,
479.791,
349.161,
409.237,
416.362,
441.355,
490.545,
436.868,

])
absorption = np.array([
0.445399,
0.341357,
0.184637,
0.154453,
0.24594,
0.126726,
0.164751,
0.169571,
0.368798,
0.222154,
0.110484,
0.0830535,
0.163056,
0.201017,
0.404513,
0.0324628,
0.373619,
])
components = fit_multimodal_gaussian(x_data=lambda_vals, y_data=absorption, n_modes=2)
for i, comp in enumerate(components):
    print(f"Mode {i+1}: height={comp['amplitude']:.3f}, mean={comp['mean']:.2f}, width={comp['width']:.2f}")
# Example usage: use plot digitizer data to get absorption parameters.
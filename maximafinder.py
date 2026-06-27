"""
Find local maxima on a sampled response surface.

The preferred workflow supports scattered samples such as Latin hypercube
sampling (LHS). A legacy regular-grid workflow is still available for complete
Cartesian grids.
"""

import json

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.stats import chi2


def build_interpolator(points, values, method="cubic"):
    """
    Build a function that estimates the response between sampled grid points.

    Parameters
    ----------
    points : list[np.ndarray]
        One sorted coordinate array per knob/dimension.
    values : np.ndarray
        Response values arranged on the grid described by points.
    method : str
        Interpolation method accepted by RegularGridInterpolator, such as
        "linear", "nearest", "slinear", "cubic", "quintic", or "pchip".

    Returns
    -------
    RegularGridInterpolator
        Callable interpolator. Pass it rows shaped like (n_samples, n_dims).
    """
    return RegularGridInterpolator(
        points,
        values,
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )


def find_grid_local_maxima(values, neighborhood_size=3):
    """
    Find candidate maxima directly on the discrete grid.

    A point is a candidate if it is equal to the largest value in a surrounding
    neighborhood. These are only coarse guesses; refine_maximum moves each one
    off the grid if the interpolated surface has a better nearby peak.
    """
    neighborhood_max = maximum_filter(values, size=neighborhood_size)
    is_local_maximum = neighborhood_max == values
    candidate_indices = np.argwhere(is_local_maximum)
    return [tuple(index) for index in candidate_indices]


def refine_maximum(interp, x0, bounds, optimizer_method="L-BFGS-B"):
    """
    Refine one grid candidate into a continuous maximum.

    scipy.optimize.minimize only minimizes, so this minimizes -response(x).
    The result is therefore the x location where response(x) is largest.
    """

    def negative_response(x):
        return -interp(x.reshape(1, -1))[0]

    result = minimize(
        negative_response,
        x0=x0,
        bounds=bounds,
        method=optimizer_method,
    )
    return result.x, -result.fun


def load_samples_from_csv(csv_path, knob_cols, value_col, chunk_size=50000):
    """
    Load a CSV as scattered samples instead of a complete grid.

    This is the right loader for LHS / Latin-hypercube data. Each row is one
    sampled design point, and the sampled points do not need to form every
    possible combination of knob values.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        samples has shape (n_samples, n_knobs). responses has shape
        (n_samples,).
    """
    print("Loading scattered samples from CSV...")

    sample_chunks = []
    response_chunks = []
    columns_to_read = knob_cols + [value_col]

    for chunk in pd.read_csv(
        csv_path,
        header=None,
        usecols=columns_to_read,
        chunksize=chunk_size,
    ):
        valid_rows = chunk.dropna()
        sample_chunks.append(valid_rows[knob_cols].to_numpy(dtype=np.float64))
        response_chunks.append(valid_rows[value_col].to_numpy(dtype=np.float64))

    if not sample_chunks:
        raise ValueError(f"No usable rows found in {csv_path}")

    samples = np.vstack(sample_chunks)
    responses = np.concatenate(response_chunks)

    print(f"Loaded {len(responses)} scattered samples with {samples.shape[1]} knobs.")
    return samples, responses


def build_scattered_model(
    samples,
    responses,
    method="idw",
    smooth=1e-6,
    rbf_function="thin_plate",
    k_neighbors=64,
    idw_power=2.0,
    max_rbf_samples=5000,
):
    """
    Build a response model from scattered samples.

    This does not require a regular grid. Coordinates are normalized to [0, 1]
    before fitting so knobs with larger numeric ranges do not dominate distance
    calculations.

    method="idw" uses local inverse-distance weighting on a KD-tree. It scales
    well for large LHS files because each evaluation only looks at k_neighbors.

    method="rbf" keeps the old radial-basis-function option for small datasets,
    but it requires an all-pairs matrix and is intentionally capped.
    """
    samples = np.asarray(samples, dtype=np.float64)
    responses = np.asarray(responses, dtype=np.float64)

    lower = samples.min(axis=0)
    upper = samples.max(axis=0)
    scale = upper - lower
    scale[scale == 0] = 1.0

    normalized_samples = (samples - lower) / scale

    if method == "idw":
        tree = cKDTree(normalized_samples)
        neighbor_count = min(k_neighbors, len(samples))

        def model(x):
            query = np.asarray(x, dtype=np.float64)
            if query.ndim == 1:
                query = query.reshape(1, -1)

            normalized_query = (query - lower) / scale
            distances, indices = tree.query(normalized_query, k=neighbor_count)

            if neighbor_count == 1:
                distances = distances.reshape(-1, 1)
                indices = indices.reshape(-1, 1)

            estimates = np.empty(len(normalized_query), dtype=np.float64)
            for row_index, (row_distances, row_indices) in enumerate(zip(distances, indices)):
                exact_matches = row_distances <= 1e-12
                if np.any(exact_matches):
                    estimates[row_index] = responses[row_indices[exact_matches][0]]
                    continue

                weights = 1.0 / np.power(row_distances + smooth, idw_power)
                estimates[row_index] = np.dot(weights, responses[row_indices]) / weights.sum()

            return estimates

        return model

    if method != "rbf":
        raise ValueError("Scattered samples support method='idw' or method='rbf'.")

    if len(samples) > max_rbf_samples:
        raise ValueError(
            "method='rbf' is too memory-hungry for this many scattered samples. "
            f"It received {len(samples)} rows, but max_rbf_samples is "
            f"{max_rbf_samples}. Use method='idw' for large LHS data."
        )

    from scipy.interpolate import Rbf

    rbf = Rbf(
        *normalized_samples.T,
        responses,
        function=rbf_function,
        smooth=smooth,
    )

    def model(x):
        query = np.asarray(x, dtype=np.float64)
        if query.ndim == 1:
            query = query.reshape(1, -1)
        normalized_query = (query - lower) / scale
        return np.asarray(rbf(*normalized_query.T))

    return model


def find_all_maxima_scattered(
    samples,
    responses,
    n_starts=50,
    method="idw",
    dedup_tol=1e-3,
    smooth=1e-6,
    rbf_function="thin_plate",
    k_neighbors=64,
    idw_power=2.0,
    optimizer_method="Powell",
):
    """
    Find maxima from scattered samples, suitable for LHS data.

    The function fits a smooth surrogate model to the sampled points, starts
    optimization from the best observed samples, and deduplicates optimized
    peaks that converge to the same location.
    """
    samples = np.asarray(samples, dtype=np.float64)
    responses = np.asarray(responses, dtype=np.float64)

    if samples.ndim != 2:
        raise ValueError("samples must have shape (n_samples, n_knobs)")
    if len(samples) != len(responses):
        raise ValueError("samples and responses must contain the same number of rows")

    model = build_scattered_model(
        samples,
        responses,
        method=method,
        smooth=smooth,
        rbf_function=rbf_function,
        k_neighbors=k_neighbors,
        idw_power=idw_power,
    )
    bounds = list(zip(samples.min(axis=0), samples.max(axis=0)))

    start_count = min(n_starts, len(responses))
    start_indices = np.argsort(responses)[-start_count:][::-1]
    print(f"Refining from the best {start_count} observed samples")

    refined_maxima = []
    for start_index in start_indices:
        start_point = samples[start_index]
        peak_location, peak_value = refine_maximum(
            model,
            start_point,
            bounds,
            optimizer_method=optimizer_method,
        )
        if np.isfinite(peak_value):
            refined_maxima.append((peak_location, peak_value))

    unique_maxima = []
    for peak_location, peak_value in sorted(refined_maxima, key=lambda item: -item[1]):
        is_new_peak = all(
            np.linalg.norm(peak_location - accepted["x"]) > dedup_tol
            for accepted in unique_maxima
        )
        if is_new_peak:
            unique_maxima.append({"x": peak_location, "f": peak_value})

    return unique_maxima, model


def find_all_maxima(
    points,
    values,
    neighborhood_size=3,
    method="cubic",
    dedup_tol=1e-3,
    smooth_sigma=1.5,
):
    """
    Find local maxima of a gridded response surface.

    Pipeline:
    1. Smooth noisy grid values with a Gaussian blur.
    2. Build an interpolator over the smoothed grid.
    3. Find coarse maximum candidates on the grid.
    4. Refine each candidate with continuous optimization.
    5. Remove duplicate peaks that converged to the same location.

    Returns
    -------
    tuple[list[dict], RegularGridInterpolator]
        maxima is a list of {"x": location, "f": response} dictionaries sorted
        from highest response to lowest. interp is the smoothed interpolator.
    """
    smoothed_values = apply_gaussian_smoothing(values, sigma=smooth_sigma)
    interp = build_interpolator(points, smoothed_values, method=method)
    bounds = [(dimension.min(), dimension.max()) for dimension in points]

    candidate_indices = find_grid_local_maxima(
        smoothed_values,
        neighborhood_size=neighborhood_size,
    )
    print(f"Found {len(candidate_indices)} grid candidates")

    refined_maxima = []
    for grid_index in candidate_indices:
        start_point = np.array(
            [points[dimension][grid_index[dimension]] for dimension in range(len(points))]
        )
        peak_location, peak_value = refine_maximum(interp, start_point, bounds)
        refined_maxima.append((peak_location, peak_value))

    unique_maxima = []
    for peak_location, peak_value in sorted(refined_maxima, key=lambda item: -item[1]):
        is_new_peak = all(
            np.linalg.norm(peak_location - accepted["x"]) > dedup_tol
            for accepted in unique_maxima
        )
        if is_new_peak:
            unique_maxima.append({"x": peak_location, "f": peak_value})

    return unique_maxima, interp


def apply_gaussian_smoothing(values, sigma=1.0):
    """
    Smooth the response grid with a Gaussian blur.

    Larger sigma values remove more small bumps/noise. Smaller values preserve
    sharper local structure.
    """
    return gaussian_filter(values, sigma=sigma)


def load_grid_from_csv_memory_efficient(csv_path, knob_cols, value_col, chunk_size=50000):
    """
    Load a CSV response table into a regular multidimensional grid.

    Use this only when the CSV contains a complete Cartesian grid. For LHS or
    other scattered samples, use load_samples_from_csv instead.

    The function reads the CSV twice to avoid building a huge temporary table:
    first it discovers the grid coordinates, then it allocates the final grid
    and fills it chunk by chunk.

    CSV columns are zero-based because pandas usecols is zero-based.
    """
    print("Pass 1: scanning CSV for grid coordinates...")

    unique_values_by_column = {column: set() for column in knob_cols}
    total_rows = 0
    columns_to_read = knob_cols + [value_col]

    for chunk in pd.read_csv(
        csv_path,
        header=None,
        usecols=columns_to_read,
        chunksize=chunk_size,
    ):
        for column in knob_cols:
            unique_values_by_column[column].update(chunk[column].unique())
        total_rows += len(chunk)

    points = [
        np.sort(np.array(list(unique_values_by_column[column])))
        for column in knob_cols
    ]
    grid_shape = tuple(len(dimension_points) for dimension_points in points)
    expected_rows = int(np.prod(grid_shape))

    print(f"Grid shape: {' x '.join(map(str, grid_shape))} = {expected_rows} points")
    print(f"CSV rows: {total_rows}")

    if total_rows != expected_rows:
        raise ValueError(
            "This CSV is not a complete Cartesian grid: "
            f"{total_rows} rows were found, but the unique knob values imply "
            f"{expected_rows} grid points. Use load_samples_from_csv and "
            "find_all_maxima_scattered for LHS/scattered data."
        )

    print(f"Allocating response grid ({expected_rows * 8 / 1e9:.2f} GB)...")
    values = np.empty(grid_shape, dtype=np.float64)

    value_to_grid_index = [
        {coordinate_value: index for index, coordinate_value in enumerate(points[dimension])}
        for dimension in range(len(knob_cols))
    ]

    print("Pass 2: filling response grid...")
    filled_rows = 0

    for chunk in pd.read_csv(
        csv_path,
        header=None,
        usecols=columns_to_read,
        chunksize=chunk_size,
    ):
        for _, row in chunk.iterrows():
            try:
                grid_index = tuple(
                    value_to_grid_index[dimension][row[knob_cols[dimension]]]
                    for dimension in range(len(knob_cols))
                )
            except KeyError:
                continue

            values[grid_index] = row[value_col]
            filled_rows += 1

        if filled_rows % (chunk_size * 5) == 0:
            print(f"  Filled {filled_rows}/{total_rows} rows...")

    print(f"Done. Loaded {filled_rows} rows into the response grid.")
    return points, values


def grid_is_complete(samples):
    """
    Return True when scattered samples form a complete Cartesian grid.

    LHS data should normally return False because it samples selected points
    from the design space rather than every knob combination.
    """
    samples = np.asarray(samples)
    unique_counts = [len(np.unique(samples[:, dimension])) for dimension in range(samples.shape[1])]
    expected_grid_rows = int(np.prod(unique_counts))
    return len(samples) == expected_grid_rows


def samples_to_regular_grid(samples, responses):
    """
    Convert scattered samples to a regular grid when they are complete.

    Raises ValueError if any grid cell is missing.
    """
    samples = np.asarray(samples, dtype=np.float64)
    responses = np.asarray(responses, dtype=np.float64)

    points = [np.unique(samples[:, dimension]) for dimension in range(samples.shape[1])]
    grid_shape = tuple(len(dimension_points) for dimension_points in points)
    values = np.full(grid_shape, np.nan, dtype=np.float64)

    value_to_grid_index = [
        {coordinate_value: index for index, coordinate_value in enumerate(points[dimension])}
        for dimension in range(samples.shape[1])
    ]

    for sample, response in zip(samples, responses):
        grid_index = tuple(
            value_to_grid_index[dimension][sample[dimension]]
            for dimension in range(samples.shape[1])
        )
        values[grid_index] = response

    if np.isnan(values).any():
        missing_count = int(np.isnan(values).sum())
        raise ValueError(f"Samples are missing {missing_count} regular-grid cells")

    return points, values


def save_maxima_to_csv(maxima, output_path):
    """Save maxima as CSV rows: knob_1, knob_2, ..., knob_n, response."""
    rows = []
    for maximum in maxima:
        rows.append(list(maximum["x"]) + [maximum["f"]])

    pd.DataFrame(rows).to_csv(output_path, index=False, header=False)
    print(f"Saved {len(maxima)} maxima to {output_path}")


def save_maxima_to_json(maxima, output_path):
    """Save maxima to a human-readable JSON file."""
    json_data = [
        {"x": maximum["x"].tolist(), "f": float(maximum["f"])}
        for maximum in maxima
    ]
    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(json_data, file, indent=2)
    print(f"Saved {len(maxima)} maxima to {output_path}")


def save_maxima_to_mat(maxima, output_path):
    """Save maxima to a MATLAB .mat file with arrays named x and f."""
    from scipy.io import savemat

    x_array = np.array([maximum["x"] for maximum in maxima])
    f_array = np.array([maximum["f"] for maximum in maxima])

    savemat(output_path, {"x": x_array, "f": f_array})
    print(f"Saved {len(maxima)} maxima to {output_path}")


def save_ellipse_to_json(
    output_path,
    center,
    cov_matrix,
    eigenvals,
    eigenvecs,
    selected_count,
    selection_summary,
    confidence=None,
    empirical_coverage=None,
):
    """Save ellipsoid parameters to a human-readable JSON file."""
    result = {
        "center": center.tolist(),
        "covariance_matrix": cov_matrix.tolist(),
        "principal_variances": eigenvals.tolist(),
        "principal_directions": eigenvecs.tolist(),
        "selected_count": int(selected_count),
        "selection": selection_summary,
    }

    if confidence is not None:
        scale = float(np.sqrt(chi2.ppf(confidence, len(center))))
        result["confidence"] = float(confidence)
        result["confidence_axis_lengths"] = (np.sqrt(eigenvals) * scale).tolist()

    if empirical_coverage is not None:
        result["empirical_coverage"] = empirical_coverage

    with open(output_path, "w", encoding="utf-8") as file:
        json.dump(result, file, indent=2)
    print(f"Saved ellipsoid parameters to {output_path}")


def fit_smooth_model(points, values, method="idw"):
    """
    Fit a smooth response model that can be evaluated away from grid points.

    This is optional and is not used by the default maximum-finding pipeline.
    The default "idw" method is scalable for scattered/LHS data. The "rbf"
    method is only practical for small datasets.
    """
    if isinstance(points, np.ndarray) and points.ndim == 2:
        return build_scattered_model(points, values, method=method)

    from scipy.interpolate import Rbf

    grid = np.meshgrid(*points, indexing="ij")
    coordinates = np.column_stack([axis_values.ravel() for axis_values in grid])
    flat_values = values.ravel()

    valid = ~np.isnan(flat_values)
    coordinates = coordinates[valid]
    flat_values = flat_values[valid]

    if method == "rbf":
        rbf = Rbf(*coordinates.T, flat_values, function="thin_plate", smooth=0.1)
        return lambda x: rbf(*x.T)

    return build_interpolator(points, values, method="cubic")


def flatten_response_data(points, values):
    """
    Return coordinate rows and response values from grid or scattered inputs.
    """
    if isinstance(points, np.ndarray) and points.ndim == 2:
        coordinates = np.asarray(points, dtype=np.float64)
        flat_values = np.asarray(values, dtype=np.float64)
    else:
        grid = np.meshgrid(*points, indexing="ij")
        coordinates = np.column_stack([axis_values.ravel() for axis_values in grid])
        flat_values = np.asarray(values, dtype=np.float64).ravel()

    valid = ~np.isnan(flat_values)
    return coordinates[valid], flat_values[valid]


def select_response_region(
    points,
    values,
    min_response=None,
    max_response_tol=None,
    top_percentile=None,
):
    """
    Select sampled points that define the high-response region.

    Use max_response_tol when the best region is a plateau, for example all
    points with f close to the maximum value 1. Use min_response when there is
    a known acceptable-response cutoff. Use top_percentile as a fallback when
    no physical threshold is known.
    """
    coordinates, flat_values = flatten_response_data(points, values)
    max_response = float(np.max(flat_values))

    if max_response_tol is not None:
        threshold = max_response - max_response_tol
        summary = {
            "mode": "near_max",
            "max_response": max_response,
            "max_response_tol": float(max_response_tol),
            "threshold": float(threshold),
        }
    elif min_response is not None:
        threshold = min_response
        summary = {
            "mode": "threshold",
            "min_response": float(min_response),
            "threshold": float(threshold),
        }
    elif top_percentile is not None:
        threshold = np.percentile(flat_values, 100 - top_percentile)
        summary = {
            "mode": "top_percentile",
            "top_percentile": float(top_percentile),
            "threshold": float(threshold),
        }
    else:
        raise ValueError(
            "Specify one of max_response_tol, min_response, or top_percentile."
        )

    selected_mask = flat_values >= threshold
    selected_points = coordinates[selected_mask]
    selected_values = flat_values[selected_mask]

    if len(selected_points) < 2:
        raise ValueError("Need at least two selected points to fit an ellipsoid.")

    print(f"Response range: [{flat_values.min():.4f}, {flat_values.max():.4f}]")
    print(f"Selection mode: {summary['mode']}")
    print(f"Selected {len(selected_points)} of {len(flat_values)} sampled points")
    print(
        "Selected response range: "
        f"[{selected_values.min():.4f}, {selected_values.max():.4f}]"
    )

    return selected_points, selected_values, summary


def fit_ellipse_to_response_region(
    points,
    values,
    min_response=None,
    max_response_tol=None,
    top_percentile=None,
):
    """
    Fit a covariance ellipsoid to a selected high-response region.

    This is the function to use when many sampled points share the maximum
    response, such as a large f = 1 plateau.
    """
    selected_points, selected_values, summary = select_response_region(
        points,
        values,
        min_response=min_response,
        max_response_tol=max_response_tol,
        top_percentile=top_percentile,
    )

    center = selected_points.mean(axis=0)
    cov_matrix = np.cov(selected_points.T)

    return center, cov_matrix, selected_points, selected_values, summary


def fit_ellipse_to_best_points(points, values, top_percentile=20):
    """
    Fit a covariance ellipsoid around the highest-response sampled points.

    This is a descriptive ellipsoid, not a formal confidence ellipsoid. It:
    1. accepts scattered rows or flattens a regular grid,
    2. keeps the top response percentile,
    3. returns the mean and covariance of those selected coordinates.
    """
    center, cov_matrix, best_points, _, _ = fit_ellipse_to_response_region(
        points,
        values,
        top_percentile=top_percentile,
    )

    return center, cov_matrix, best_points


def ellipse_empirical_coverage(points, center, cov_matrix, confidence):
    """
    Measure what fraction of actual sampled points fall inside an ellipsoid.

    This is different from the theoretical Gaussian coverage. It directly
    counts the selected data points whose squared Mahalanobis distance is below
    the chi-square threshold for the requested confidence level.
    """
    points = np.asarray(points, dtype=np.float64)
    center = np.asarray(center, dtype=np.float64)

    threshold = float(chi2.ppf(confidence, len(center)))
    inv_cov = np.linalg.pinv(cov_matrix)
    deltas = points - center
    mahalanobis_squared = np.einsum("ij,jk,ik->i", deltas, inv_cov, deltas)
    inside = mahalanobis_squared <= threshold

    inside_count = int(np.count_nonzero(inside))
    total_count = int(len(points))
    fraction = inside_count / total_count

    print(
        f"Empirical coverage of selected points inside {confidence:.1%} ellipsoid: "
        f"{inside_count}/{total_count} = {100 * fraction:.2f}%"
    )

    return {
        "confidence": float(confidence),
        "chi_square_threshold": threshold,
        "inside_count": inside_count,
        "total_count": total_count,
        "fraction": float(fraction),
        "percent": float(100 * fraction),
    }


def describe_ellipse(center, cov_matrix, confidence=None):
    """
    Print the center, covariance, and principal axes of an ellipsoid.

    Eigenvectors are the ellipsoid directions. The square root of each
    eigenvalue is the one-standard-deviation length along that direction.
    If confidence is provided, also print chi-square-scaled axis lengths.
    """
    print("\n=== HYPERELLIPSE FROM HIGH-RESPONSE POINTS ===")
    print(f"Center (mean): {np.round(center, 4)}")
    print("\nCovariance matrix:")
    print(np.round(cov_matrix, 6))

    eigenvals, eigenvecs = np.linalg.eigh(cov_matrix)
    print(f"\nPrincipal variances: {np.round(eigenvals, 6)}")
    one_sigma_axes = np.sqrt(eigenvals)
    print(f"Principal std devs: {np.round(one_sigma_axes, 4)}")

    if confidence is not None:
        scale = np.sqrt(chi2.ppf(confidence, len(center)))
        print(
            f"{confidence:.1%} ellipsoid axis lengths: "
            f"{np.round(one_sigma_axes * scale, 4)}"
        )

    print("\nPrincipal directions:")
    for index, eigenvector in enumerate(eigenvecs.T):
        print(f"  Direction {index}: {np.round(eigenvector, 4)}")

    return eigenvals, eigenvecs


def estimate_ellipse_integral(center, cov_matrix, model_func, num_samples=10000):
    """
    Estimate the response integral near the ellipsoid with Monte Carlo samples.

    Samples are drawn from the Gaussian described by center and cov_matrix.
    The returned value is a rough estimate of mean response times the Gaussian
    one-sigma volume scale.
    """
    samples = np.random.multivariate_normal(center, cov_matrix, num_samples)
    evaluations = model_func(samples)

    det_cov = np.linalg.det(cov_matrix)
    volume = np.sqrt((2 * np.pi) ** len(center) * det_cov)
    integral_estimate = np.mean(evaluations) * volume

    print(f"\nMonte Carlo integral estimate: {integral_estimate:.6f}")
    print(f"Mean response near ellipsoid: {np.mean(evaluations):.6f}")
    print(f"Effective one-sigma volume: {volume:.6f}")

    return integral_estimate


if __name__ == "__main__":
    # ---- User settings -------------------------------------------------
    # CSV to analyze. The file is expected to have no header row.
    csv_path = "IDAa_res9U_N.csv"

    # Zero-based CSV columns that define the grid coordinates.
    knob_cols = [6, 7, 8, 9, 10, 11]

    # Zero-based CSV column containing the response to maximize.
    value_col = 12

    # Output file for the detected maxima.
    output_csv = "IDAa_cons9U_N_maxima_results.csv"
    output_ellipse_json = "IDAa_cons9U_N_ellipsoid.json"

    # Data mode:
    # "scattered" is the right choice for LHS data and does not require a full grid.
    # "grid" requires a complete Cartesian grid.
    # "auto" uses the grid path only when the samples form a complete grid.
    data_mode = "auto"

    # Shared controls.
    chunk_size = 50000
    dedup_tol = 1e-3
    run_maxima_search = False

    # Scattered/LHS surrogate controls.
    n_starts = 30
    scattered_method = "idw"
    k_neighbors = 64
    idw_power = 2.0
    optimizer_method = "Powell"
    rbf_function = "thin_plate"
    rbf_smooth = 1e-6

    # Regular-grid-only controls.
    neighborhood_size = 5
    interpolation_method = "pchip"
    smooth_sigma = 1.5

    # Ellipsoid summary controls. Works for both scattered and grid data.
    # "near_max": use all points within max_response_tol of the observed maximum.
    # "threshold": use all points with response >= ellipse_min_response.
    # "top_percentile": use the highest top_percentile percent of responses.
    ellipse_selection = "near_max"
    max_response_tol = 1e-12
    ellipse_min_response = 1.0
    top_percentile = 70
    confidence = 0.95
    # ------------------------------------------------------------------

    samples, responses = load_samples_from_csv(
        csv_path,
        knob_cols,
        value_col,
        chunk_size=chunk_size,
    )

    use_grid_path = data_mode == "grid" or (
        data_mode == "auto" and grid_is_complete(samples)
    )

    if use_grid_path:
        print("Using complete-grid maximum finder.")
        points, response_grid = samples_to_regular_grid(samples, responses)
        ellipse_points = points
        ellipse_values = response_grid
    else:
        print("Using scattered/LHS maximum finder.")
        ellipse_points = samples
        ellipse_values = responses

    if run_maxima_search:
        if use_grid_path:
            maxima, model = find_all_maxima(
                points,
                response_grid,
                neighborhood_size=neighborhood_size,
                method=interpolation_method,
                dedup_tol=dedup_tol,
                smooth_sigma=smooth_sigma,
            )
        else:
            maxima, model = find_all_maxima_scattered(
                samples,
                responses,
                n_starts=n_starts,
                method=scattered_method,
                dedup_tol=dedup_tol,
                smooth=rbf_smooth,
                rbf_function=rbf_function,
                k_neighbors=k_neighbors,
                idw_power=idw_power,
                optimizer_method=optimizer_method,
            )

        print(f"\nFound {len(maxima)} local maxima")
        save_maxima_to_csv(maxima, output_csv)
    else:
        print("Skipping local-maxima search; fitting ellipsoid to response region.")

    if ellipse_selection == "near_max":
        ellipse_kwargs = {"max_response_tol": max_response_tol}
    elif ellipse_selection == "threshold":
        ellipse_kwargs = {"min_response": ellipse_min_response}
    elif ellipse_selection == "top_percentile":
        ellipse_kwargs = {"top_percentile": top_percentile}
    else:
        raise ValueError(
            "ellipse_selection must be 'near_max', 'threshold', or 'top_percentile'"
        )

    center, cov_matrix, best_points, best_values, selection_summary = (
        fit_ellipse_to_response_region(
            ellipse_points,
            ellipse_values,
            **ellipse_kwargs,
        )
    )
    eigenvals, eigenvecs = describe_ellipse(
        center,
        cov_matrix,
        confidence=confidence,
    )
    empirical_coverage = ellipse_empirical_coverage(
        best_points,
        center,
        cov_matrix,
        confidence,
    )
    save_ellipse_to_json(
        output_ellipse_json,
        center,
        cov_matrix,
        eigenvals,
        eigenvecs,
        len(best_points),
        selection_summary,
        confidence=confidence,
        empirical_coverage=empirical_coverage,
    )

    # Backward-compatible percentile ellipsoid call:
    # center, cov_matrix, best_points = fit_ellipse_to_best_points(
    #     ellipse_points,
    #     ellipse_values,
    #     top_percentile=top_percentile,
    # )
    # eigenvals, eigenvecs = describe_ellipse(center, cov_matrix)

    # Optional integral estimate:
    # model = fit_smooth_model(ellipse_points, ellipse_values, method="idw")
    # estimate_ellipse_integral(center, cov_matrix, model, num_samples=10000)

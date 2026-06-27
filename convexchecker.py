import numpy as np
import pandas as pd


def load_regular_grid(csv_path, grid_cols, value_col, delimiter=","):
    df = pd.read_csv(csv_path, delimiter=delimiter, header=None)

    points = [np.sort(df[c].unique()) for c in grid_cols]
    shape = tuple(len(p) for p in points)
    expected_rows = int(np.prod(shape))

    if len(df) != expected_rows:
        raise ValueError(
            f"CSV is not a full regular grid: expected {expected_rows} rows, got {len(df)}"
        )

    df_sorted = df.sort_values(grid_cols)
    values = df_sorted[value_col].to_numpy().reshape(shape)
    return points, values


def check_uniform_spacing(points, tol=1e-8):
    for i, axis in enumerate(points):
        diffs = np.diff(axis)
        if not np.allclose(diffs, diffs[0], atol=tol, rtol=0):
            raise ValueError(
                f"Axis {i} is not uniformly spaced; convexity check requires a regular grid"
            )
    return [np.diff(axis)[0] for axis in points]


def _shift_index(idx, axis, delta):
    idx = list(idx)
    idx[axis] += delta
    return tuple(idx)


def is_convex_on_grid(points, values, tol=1e-8):
    n = len(points)
    shape = values.shape

    if any(s < 3 for s in shape):
        raise ValueError("Each grid axis must have at least 3 points for Hessian estimation")

    spacings = check_uniform_spacing(points, tol=tol)

    failures = []
    for center in np.ndindex(*shape):
        if any(i in (0, shape[ax] - 1) for ax, i in enumerate(center)):
            continue

        H = np.zeros((n, n), dtype=float)

        for i in range(n):
            h = spacings[i]
            f_plus = values[_shift_index(center, i, +1)]
            f_minus = values[_shift_index(center, i, -1)]
            f0 = values[center]
            H[i, i] = (f_plus + f_minus - 2 * f0) / (h * h)

        for i in range(n):
            for j in range(i + 1, n):
                hi = spacings[i]
                hj = spacings[j]
                f_pp = values[_shift_index(_shift_index(center, i, +1), j, +1)]
                f_pm = values[_shift_index(_shift_index(center, i, +1), j, -1)]
                f_mp = values[_shift_index(_shift_index(center, i, -1), j, +1)]
                f_mm = values[_shift_index(_shift_index(center, i, -1), j, -1)]
                H[i, j] = H[j, i] = (f_pp - f_pm - f_mp + f_mm) / (4 * hi * hj)

        eigs = np.linalg.eigvalsh(H)
        if eigs.min() < -tol:
            failures.append(
                {"point": center, "min_eig": float(eigs.min()), "hessian": H}
            )

    return len(failures) == 0, failures


def is_concave_on_grid(points, values, tol=1e-8):
    return is_convex_on_grid(points, -values, tol=tol)


def main():
    csv_path = "IDAa_res9U_P.csv"
    grid_cols = [6, 7, 8, 9, 10, 11]  # adapt to your dimensionality
    value_col = 12               # the function column

    points, values = load_regular_grid(csv_path, grid_cols, value_col)
    concave, failures = is_concave_on_grid(points, values, tol=1e-6)

    if concave:
        print("The sampled function appears concave on the regular grid.")
    else:
        print("The function is not concave on the sampled grid.")
        print(f"Found {len(failures)} violating interior points.")
        for fail in failures[:10]:
            print(f"  index={fail['point']} min_eig={fail['min_eig']:.3e}")


if __name__ == "__main__":
    main()
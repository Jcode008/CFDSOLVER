import pandas as pd

# Load metrics
metrics = pd.read_csv('build/metrics.csv')

# Check point (74, 1)
i, j = 74, 1
nxi = 200
idx = j * nxi + i

print(f"Metrics at (i={i}, j={j}):")
print(f"  J     = {metrics.loc[idx, 'J']}")
print(f"  xi_x  = {metrics.loc[idx, 'xi_x']}")
print(f"  xi_y  = {metrics.loc[idx, 'xi_y']}")
print(f"  eta_x = {metrics.loc[idx, 'eta_x']}")
print(f"  eta_y = {metrics.loc[idx, 'eta_y']}")

# Check neighbors
print(f"\nNeighbor Jacobians:")
print(f"  J[{i-1},{j}] = {metrics.loc[j*nxi+(i-1), 'J']}")
print(f"  J[{i+1},{j}] = {metrics.loc[j*nxi+(i+1), 'J']}")
print(f"  J[{i},{j-1}] = {metrics.loc[(j-1)*nxi+i, 'J']}")
print(f"  J[{i},{j+1}] = {metrics.loc[(j+1)*nxi+i, 'J']}")

# Check for extremely small/large values
print(f"\nMetric magnitude checks:")
print(f"  alpha = xi_x^2 + xi_y^2 = {metrics.loc[idx, 'xi_x']**2 + metrics.loc[idx, 'xi_y']**2}")
print(f"  gamma = eta_x^2 + eta_y^2 = {metrics.loc[idx, 'eta_x']**2 + metrics.loc[idx, 'eta_y']**2}")

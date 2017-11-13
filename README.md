# ez-shocks

Adapted from Jurado, Ludvigson, Ng (2015) using code supplied by Serena Ng at:

- http://www.columbia.edu/~sn2294/research.html

## Generate macro uncertainty

1. Clean and suitably transform data.

```
import_data.m
```

2. Estimate number of latent factors in data and extract via PCA. Forecast macro variables and save prediction errors.

```
factors_forecast.m
```

3. Estimate latent stochastic volatility in prediction errors

```
svydraws.R
svfdraws.R
```

4. Compute uncertainty in variables and build aggregate index

```
uncertainty.m
```

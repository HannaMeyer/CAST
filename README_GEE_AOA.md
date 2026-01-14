# AOA for Google Earth Engine - Quick Start Guide

## Overview

This implementation brings the **Area of Applicability (AOA)** methodology from the CAST R package to Google Earth Engine Python API. It's specifically designed to handle class-imbalanced classification problems like mowing/noMowing detection.

## Key Features

✅ **Class-specific AOA** - Calculate AOA for specific classes (e.g., only mowing areas)
✅ **Handles class imbalance** - Don't let abundant noMowing samples dominate your AOA
✅ **GEE-native** - Works with ee.Image and ee.FeatureCollection
✅ **Spatial CV support** - Respects spatial cross-validation folds
✅ **Variable weighting** - Uses variable importance from your model
✅ **Large-scale compatible** - Optimized for continent-scale analysis

## The Problem This Solves

**Your situation:**
- You have training data with **many noMowing samples** and **few mowing samples**
- You trained a model with spatial CV (one region out) - F1 scores look good
- You tested on Slovakia & Swiss - F1 scores even better!
- Now you want to predict on **entire Alps & Carpathians** with uncertainty estimates

**The issue:**
If you calculate standard AOA with imbalanced data, it will show "inside AOA" wherever the prediction area is similar to **any training sample** - mostly noMowing areas. This doesn't help you understand where your model can reliably detect **mowing**.

**The solution:**
Calculate **class-specific AOA using only mowing samples as reference**. This shows where the prediction area is similar to known mowing areas, giving you true uncertainty estimates for mowing detection.

## Installation

1. Copy files to your project:
   ```bash
   # Core implementation
   gee_aoa.py

   # Example for mowing classification
   example_mowing_aoa.py

   # This guide
   README_GEE_AOA.md
   ```

2. Initialize Earth Engine:
   ```python
   import ee
   ee.Initialize()
   ```

## Quick Start: 3-Minute Example

```python
import ee
from gee_aoa import calculate_aoa_for_class, AOA, create_training_fc

# Initialize
ee.Initialize()

# 1. Prepare your training data (replace with actual data)
training_data = [
    {'longitude': 10.5, 'latitude': 47.5, 'NDVI': 0.8, 'elevation': 1500,
     'class': 'mowing', 'region': 'region_1'},
    {'longitude': 10.6, 'latitude': 47.6, 'NDVI': 0.3, 'elevation': 800,
     'class': 'noMowing', 'region': 'region_1'},
    # ... many more samples, mostly noMowing
]

# 2. Create feature collection
train_fc = create_training_fc(training_data, ['NDVI', 'elevation'])

# 3. Calculate AOA for MOWING ONLY (this is the key!)
mowing_aoa_params = calculate_aoa_for_class(
    train_fc=train_fc,
    predictor_names=['NDVI', 'elevation'],
    class_column='class',
    target_class='mowing',  # Only use mowing samples!
    weights={'NDVI': 0.7, 'elevation': 0.3}  # From your model
)

# 4. Load your prediction image
prediction_image = ee.Image('users/your_username/alps_carpathians_composite')

# 5. Calculate AOA
aoa_calculator = AOA(mowing_aoa_params)
result = aoa_calculator.calculate_di_optimized(prediction_image, scale=30)

# 6. Visualize
print(f"DI Threshold: {mowing_aoa_params.threshold:.4f}")

# In GEE Code Editor or geemap:
# Map.addLayer(result.select('DI'), {'min': 0, 'max': 2, 'palette': ['green', 'yellow', 'red']}, 'DI')
# Map.addLayer(result.select('AOA'), {'min': 0, 'max': 1, 'palette': ['red', 'green']}, 'AOA')
```

## Full Workflow for Mowing Classification

See `example_mowing_aoa.py` for a complete, documented workflow including:

1. ✅ Loading training data with spatial CV folds
2. ✅ Calculating standard vs. class-specific AOA
3. ✅ Applying AOA to Alps & Carpathians
4. ✅ Comparing results and visualization
5. ✅ Exporting to GEE assets
6. ✅ Validating test data (Slovakia/Swiss) against AOA
7. ✅ Interpretation guide for decision-making

Run it:
```python
python example_mowing_aoa.py
```

## Key Concepts

### Dissimilarity Index (DI)

The DI measures how similar a new location is to your training data:

- **Low DI (< threshold)**: Similar to training mowing areas → Reliable predictions
- **High DI (> threshold)**: Dissimilar to training mowing areas → Uncertain predictions

Formula: `DI = min_distance_to_training / avg_training_distance`

### Area of Applicability (AOA)

Binary classification based on DI:

- **AOA = 1 (Inside)**: DI ≤ threshold → Model can reliably predict here
- **AOA = 0 (Outside)**: DI > threshold → Predictions are uncertain

Threshold is calculated as: `Q3 + 1.5 × IQR` of training data DI values

### Why Class-Specific Matters

**Standard AOA (using all training data):**
```
Training: 1000 noMowing + 50 mowing samples
Reference: All 1050 samples
Result: AOA = 1 wherever similar to ANY training sample
Problem: Dominated by abundant noMowing areas
```

**Class-Specific AOA (using only mowing):**
```
Training: 1000 noMowing + 50 mowing samples
Reference: Only 50 mowing samples
Result: AOA = 1 wherever similar to MOWING areas
Benefit: True uncertainty for mowing detection
```

## Files Description

| File | Purpose |
|------|---------|
| `gee_aoa.py` | Core implementation (TrainDI, AOA classes) |
| `example_mowing_aoa.py` | Complete workflow for mowing classification |
| `README_GEE_AOA.md` | This guide |

## API Reference

### Class: `TrainDI`

Calculates training statistics for AOA.

**Constructor:**
```python
TrainDI(train, predictor_names, weights=None, cv_folds=None, method="L2", use_cv=False)
```

**Parameters:**
- `train`: ee.FeatureCollection or list of dicts with training data
- `predictor_names`: List of predictor variable names
- `weights`: Dict of variable importance (default: all 1.0)
- `cv_folds`: List of fold assignments for spatial CV (optional)
- `method`: Distance metric ("L2" for Euclidean)
- `use_cv`: Use CV folds for threshold calculation

**Key Attributes:**
- `threshold`: DI threshold for AOA determination
- `train_dist_avrgmean`: Normalization factor
- `scale_mean`, `scale_std`: Scaling parameters
- `weights`: Variable importance weights

**Methods:**
- `to_dict()`: Export parameters as dictionary
- `from_dict(params)`: Load from dictionary

### Class: `AOA`

Applies AOA to prediction data.

**Constructor:**
```python
AOA(train_di)
```

**Parameters:**
- `train_di`: TrainDI object with training statistics

**Methods:**

#### `calculate_di_optimized(image, geometry=None, scale=30, tile_scale=1)`

Calculate DI and AOA for an image.

**Parameters:**
- `image`: ee.Image with predictor bands
- `geometry`: Region of interest (optional)
- `scale`: Resolution in meters (default: 30)
- `tile_scale`: Tiling factor for large areas (1-16)

**Returns:**
- ee.Image with 'DI' and 'AOA' bands

### Function: `calculate_aoa_for_class()`

Calculate class-specific AOA.

```python
calculate_aoa_for_class(
    train_fc,           # Training feature collection
    predictor_names,    # List of predictors
    class_column,       # Column with class labels
    target_class,       # Class to use (e.g., "mowing")
    weights=None,       # Variable importance
    cv_folds=None       # CV fold assignments
)
```

**Returns:** TrainDI object for target class only

### Function: `create_training_fc()`

Convert training data to ee.FeatureCollection.

```python
create_training_fc(
    data,              # List of dicts
    predictor_names,   # List of predictors
    lat_col='latitude',
    lon_col='longitude'
)
```

**Returns:** ee.FeatureCollection

## Interpreting Results

### DI Values

| DI Range | Interpretation | Action |
|----------|----------------|--------|
| 0.0 - 0.5 | Very similar to training | High confidence predictions |
| 0.5 - threshold | Similar enough | Moderate confidence |
| threshold - 2.0 | Dissimilar | Low confidence, use with caution |
| > 2.0 | Very dissimilar | Unreliable, consider as "unknown" |

### AOA for Decision Making

**Scenario 1: Model predicts "mowing", AOA = 1**
- ✅ High confidence
- Area is similar to training mowing sites
- Trust this prediction

**Scenario 2: Model predicts "mowing", AOA = 0**
- ⚠️ Low confidence
- Area is different from training mowing sites
- Might be false positive - validate in field

**Scenario 3: Model predicts "noMowing", AOA = 0**
- 🤔 Uncertain
- Can't rule out mowing in novel conditions
- Consider additional data collection

## Performance Considerations

### For Small Training Sets (< 1000 samples)

```python
# Fast and simple
train_di = TrainDI(train=train_fc, predictor_names=predictors, weights=weights)
aoa_calc = AOA(train_di)
result = aoa_calc.calculate_di_optimized(image, scale=30)
```

### For Large Training Sets (> 1000 samples)

```python
# Use tiling for large areas
result = aoa_calc.calculate_di_optimized(
    image=prediction_image,
    geometry=region,
    scale=30,
    tile_scale=4  # Process in tiles
)
```

### For Very Large Prediction Areas

```python
# Export in chunks
regions = [alps, carpathians_west, carpathians_east]
for i, region in enumerate(regions):
    result = aoa_calc.calculate_di_optimized(image, region, scale=30)
    export_di_to_asset(result, f'users/username/aoa_chunk_{i}', region)
```

## Common Issues & Solutions

### Issue: "DI values all very high"
**Cause:** Training and prediction data might be from different sensors or processing
**Solution:** Ensure same preprocessing (scaling, atmospheric correction) for both

### Issue: "AOA covers < 10% of prediction area"
**Cause:** Training data doesn't represent prediction area diversity
**Solution:**
1. Check if training covers full predictor range
2. Consider collecting more training data
3. Review if predictors are appropriate

### Issue: "Memory error with large training set"
**Cause:** Distance calculation for many training points
**Solution:**
1. Subsample training data (stratified by class)
2. Use tile_scale parameter
3. Process prediction area in chunks

### Issue: "Results differ from R CAST"
**Cause:** Small numerical differences in distance calculations
**Solution:** Differences should be < 1% - verify same:
- Scaling method (check mean/std)
- Distance metric (L2)
- Threshold calculation (Q3 + 1.5*IQR)

## Advanced: Combining with Model Predictions

```python
# Load your classification model predictions
classification = ee.Image('users/username/mowing_classification')

# Load AOA
aoa_result = ee.Image('users/username/aoa_mowing')

# Create confidence-weighted prediction
# Only show mowing predictions where AOA = 1
confident_mowing = classification.updateMask(aoa_result.select('AOA'))

# Or create 3-class output: mowing (confident), mowing (uncertain), noMowing
def create_confidence_classes(class_img, aoa_img):
    # 0: noMowing
    # 1: mowing (inside AOA)
    # 2: mowing (outside AOA - uncertain)

    is_mowing = class_img.eq(1)  # Assuming 1 = mowing
    inside_aoa = aoa_img.select('AOA').eq(1)

    # Create confidence classes
    result = class_img.where(
        is_mowing.And(inside_aoa), 1
    ).where(
        is_mowing.And(inside_aoa.Not()), 2
    )

    return result

confidence_map = create_confidence_classes(classification, aoa_result)
```

## References

1. **Original AOA Paper:**
   Meyer, H., Pebesma, E. (2021): Predicting into unknown space? Estimating the area of applicability of spatial prediction models. *Methods in Ecology and Evolution* 12: 1620-1633. https://doi.org/10.1111/2041-210X.13650

2. **CAST R Package:**
   https://github.com/HannaMeyer/CAST

3. **Google Earth Engine:**
   https://earthengine.google.com/

## Citation

If you use this implementation in your research, please cite both the original AOA paper (Meyer & Pebesma 2021) and the CAST R package.

## License

GPL (>= 2) - matching the CAST R package license

## Support

For issues specific to:
- **GEE implementation**: Check `example_mowing_aoa.py` for detailed workflow
- **AOA methodology**: See CAST R package documentation and Meyer & Pebesma (2021)
- **GEE API**: See Earth Engine documentation

## Example Output Interpretation

After running the analysis, you'll see output like:

```
Calculating mowing-specific AOA (using only mowing samples)...
  - Threshold: 0.8543
  - Training DI range: 0.2341 - 0.8543
  - Number of reference samples: 87

Applying Mowing-specific AOA to Alps and Carpathians...
  DI Statistics:
    - 5th percentile: 0.3215
    - Median: 0.6789
    - 95th percentile: 1.2456
    - Area inside AOA: 67.3%
```

**Interpretation:**
- 67.3% of Alps/Carpathians is similar enough to training mowing areas
- Median DI of 0.68 is below threshold (0.85) - good overall coverage
- 95th percentile of 1.25 shows some areas are quite dissimilar
- These dissimilar areas (32.7%) need careful interpretation or more training data

## Next Steps

1. ✅ Run `example_mowing_aoa.py` with your data
2. ✅ Compare standard vs. class-specific AOA
3. ✅ Validate that Slovakia/Swiss test areas are inside AOA
4. ✅ Export results for use in your final classification
5. ✅ Use DI values to prioritize field validation locations
6. ✅ Consider collecting more training data in high-DI areas

---

**Happy AOA Mapping!** 🗺️🌍

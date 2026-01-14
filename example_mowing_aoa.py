"""
Example: AOA for Mowing Classification in Alps and Carpathians
===============================================================

This script demonstrates how to apply Area of Applicability (AOA) analysis
to a mowing/noMowing classification problem with class imbalance.

Workflow:
1. Training data: Spatial cross-validation with one region out
2. Independent test: Slovakia and Swiss data
3. Prediction area: Entire Alps and Carpathians
4. Problem: Class imbalance (many noMowing, few mowing samples)
5. Solution: Class-specific AOA for mowing areas

Author: Based on CAST R package methodology
"""

import ee
import numpy as np
from gee_aoa import TrainDI, AOA, calculate_aoa_for_class, create_training_fc, export_di_to_asset

# Initialize Earth Engine
ee.Initialize()


# =============================================================================
# STEP 1: Prepare Training Data
# =============================================================================

def load_training_data():
    """
    Load or create training data with mowing/noMowing labels.

    In practice, you would load this from your actual training dataset.
    This example shows the expected format.
    """
    # Example structure - replace with your actual data loading
    training_data = [
        # Mowing samples (fewer)
        {'longitude': 10.5, 'latitude': 47.5, 'NDVI': 0.82, 'EVI': 0.75,
         'elevation': 1500, 'slope': 15, 'aspect': 180,
         'class': 'mowing', 'region': 'region_1'},
        {'longitude': 10.6, 'latitude': 47.6, 'NDVI': 0.78, 'EVI': 0.72,
         'elevation': 1450, 'slope': 12, 'aspect': 170,
         'class': 'mowing', 'region': 'region_1'},
        # ... more mowing samples

        # NoMowing samples (many more)
        {'longitude': 10.7, 'latitude': 47.7, 'NDVI': 0.25, 'EVI': 0.20,
         'elevation': 800, 'slope': 5, 'aspect': 90,
         'class': 'noMowing', 'region': 'region_1'},
        {'longitude': 10.8, 'latitude': 47.8, 'NDVI': 0.15, 'EVI': 0.12,
         'elevation': 750, 'slope': 3, 'aspect': 85,
         'class': 'noMowing', 'region': 'region_2'},
        # ... many more noMowing samples
    ]

    return training_data


def load_predictor_names():
    """
    Define predictor variables used in the model.
    These should match the bands in your prediction image.
    """
    return [
        'NDVI',
        'EVI',
        'elevation',
        'slope',
        'aspect',
        # Add more predictors as needed
        # 'temperature',
        # 'precipitation',
        # 'soil_moisture',
        # etc.
    ]


def load_variable_importance():
    """
    Load variable importance weights from your trained model.

    In practice, extract these from your Random Forest, XGBoost, or other model.
    Higher values = more important variables.
    """
    # Example weights - replace with actual importance from your model
    weights = {
        'NDVI': 0.35,
        'EVI': 0.28,
        'elevation': 0.15,
        'slope': 0.12,
        'aspect': 0.10,
    }

    return weights


def create_cv_folds(training_data):
    """
    Create spatial cross-validation folds based on regions.

    This matches your "one region out" cross-validation approach.
    """
    # Extract region IDs
    regions = [sample['region'] for sample in training_data]
    unique_regions = sorted(set(regions))

    # Map regions to fold numbers
    region_to_fold = {region: i for i, region in enumerate(unique_regions)}
    cv_folds = [region_to_fold[region] for region in regions]

    return cv_folds


# =============================================================================
# STEP 2: Calculate Standard AOA (for comparison)
# =============================================================================

def calculate_standard_aoa(train_fc, predictor_names, weights, cv_folds):
    """
    Calculate standard AOA using ALL training data.

    Problem: With class imbalance, this will be dominated by noMowing samples.
    """
    print("Calculating standard AOA (using all training data)...")

    train_di_all = TrainDI(
        train=train_fc,
        predictor_names=predictor_names,
        weights=weights,
        cv_folds=cv_folds,
        use_cv=True
    )

    print(f"  - Threshold: {train_di_all.threshold:.4f}")
    print(f"  - Training DI range: {min(train_di_all.train_di):.4f} - {max(train_di_all.train_di):.4f}")

    return train_di_all


# =============================================================================
# STEP 3: Calculate Class-Specific AOA for Mowing
# =============================================================================

def calculate_mowing_specific_aoa(train_fc, predictor_names, weights, cv_folds):
    """
    Calculate AOA using ONLY mowing samples as reference.

    Solution: This ensures AOA represents similarity to mowing areas specifically,
    not to the overall (noMowing-dominated) training distribution.
    """
    print("\nCalculating mowing-specific AOA (using only mowing samples)...")

    train_di_mowing = calculate_aoa_for_class(
        train_fc=train_fc,
        predictor_names=predictor_names,
        class_column='class',
        target_class='mowing',  # ONLY mowing samples used
        weights=weights,
        cv_folds=cv_folds
    )

    print(f"  - Threshold: {train_di_mowing.threshold:.4f}")
    print(f"  - Training DI range: {min(train_di_mowing.train_di):.4f} - {max(train_di_mowing.train_di):.4f}")
    print(f"  - Number of reference samples: {len(train_di_mowing.train_scaled)}")

    return train_di_mowing


# =============================================================================
# STEP 4: Load Prediction Area (Alps and Carpathians)
# =============================================================================

def load_prediction_image():
    """
    Load the image composite for Alps and Carpathians.

    Replace this with your actual image loading logic.
    The image should have bands matching your predictor_names.
    """
    # Option 1: Load from existing asset
    # prediction_image = ee.Image('users/your_username/alps_carpathians_composite')

    # Option 2: Create composite from Sentinel-2 or other data
    # Example for Sentinel-2:
    region = ee.Geometry.Rectangle([5, 43, 27, 49])  # Alps and Carpathians extent

    # Load Sentinel-2 data (example for growing season)
    s2 = ee.ImageCollection('COPERNICUS/S2_SR') \
        .filterBounds(region) \
        .filterDate('2023-05-01', '2023-09-30') \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))

    # Calculate vegetation indices
    def add_indices(image):
        ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI')
        evi = image.expression(
            '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
            {
                'NIR': image.select('B8'),
                'RED': image.select('B4'),
                'BLUE': image.select('B2')
            }
        ).rename('EVI')
        return image.addBands([ndvi, evi])

    s2_with_indices = s2.map(add_indices)

    # Create median composite
    composite = s2_with_indices.median()

    # Add topography from SRTM
    dem = ee.Image('USGS/SRTMGL1_003')
    terrain = ee.Terrain.products(dem)

    # Combine all predictors
    prediction_image = composite.select(['NDVI', 'EVI']) \
        .addBands(dem.rename('elevation')) \
        .addBands(terrain.select('slope')) \
        .addBands(terrain.select('aspect'))

    return prediction_image, region


# =============================================================================
# STEP 5: Apply AOA to Prediction Area
# =============================================================================

def apply_aoa_to_alps_carpathians(train_di, prediction_image, region, output_name):
    """
    Calculate DI and AOA for the entire Alps and Carpathians region.
    """
    print(f"\nApplying {output_name} to Alps and Carpathians...")

    aoa_calculator = AOA(train_di)

    # Calculate DI and AOA
    result = aoa_calculator.calculate_di_optimized(
        image=prediction_image,
        geometry=region,
        scale=30,  # 30m resolution
        tile_scale=4  # Use tiling for large areas
    )

    # Calculate statistics
    stats = result.select('DI').reduceRegion(
        reducer=ee.Reducer.percentile([5, 25, 50, 75, 95]),
        geometry=region,
        scale=1000,  # Coarser scale for statistics
        maxPixels=1e9
    ).getInfo()

    print(f"  DI Statistics:")
    print(f"    - 5th percentile: {stats.get('DI_p5', 0):.4f}")
    print(f"    - Median: {stats.get('DI_p50', 0):.4f}")
    print(f"    - 95th percentile: {stats.get('DI_p95', 0):.4f}")

    # Calculate AOA coverage
    aoa_stats = result.select('AOA').reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=region,
        scale=1000,
        maxPixels=1e9
    ).getInfo()

    aoa_percentage = aoa_stats.get('AOA', 0) * 100
    print(f"    - Area inside AOA: {aoa_percentage:.1f}%")

    return result


# =============================================================================
# STEP 6: Compare Results and Visualize
# =============================================================================

def compare_and_visualize(result_standard, result_mowing, region, prediction_image):
    """
    Compare standard AOA vs. mowing-specific AOA and create visualizations.
    """
    print("\n" + "="*70)
    print("COMPARISON: Standard AOA vs. Mowing-Specific AOA")
    print("="*70)

    # Calculate difference in AOA
    aoa_diff = result_mowing.select('AOA').subtract(result_standard.select('AOA'))

    # Areas that are:
    # +1: Inside mowing AOA but outside standard AOA (good for mowing detection)
    # 0: Same in both
    # -1: Inside standard AOA but outside mowing AOA (reliable in general but not for mowing)

    # Visualization parameters
    di_vis_standard = {
        'bands': ['DI'],
        'min': 0,
        'max': 2,
        'palette': ['blue', 'cyan', 'yellow', 'orange', 'red']
    }

    di_vis_mowing = {
        'bands': ['DI'],
        'min': 0,
        'max': 2,
        'palette': ['darkgreen', 'green', 'yellow', 'orange', 'red']
    }

    aoa_vis = {
        'bands': ['AOA'],
        'min': 0,
        'max': 1,
        'palette': ['red', 'green']
    }

    diff_vis = {
        'min': -1,
        'max': 1,
        'palette': ['red', 'white', 'green']
    }

    # Print visualization code (for use in Code Editor or geemap)
    print("\nVisualization code (use in GEE Code Editor or geemap):")
    print("=" * 70)
    print("""
# Add to map
Map.centerObject(region, 7)

# Standard AOA
Map.addLayer(result_standard.select('DI'), di_vis_standard, 'DI (Standard)')
Map.addLayer(result_standard.select('AOA'), aoa_vis, 'AOA (Standard)')

# Mowing-specific AOA
Map.addLayer(result_mowing.select('DI'), di_vis_mowing, 'DI (Mowing-specific)')
Map.addLayer(result_mowing.select('AOA'), aoa_vis, 'AOA (Mowing-specific)')

# Difference
Map.addLayer(aoa_diff, diff_vis, 'AOA Difference (Mowing - Standard)')

# Add prediction image for context
Map.addLayer(prediction_image.select('NDVI'),
             {'min': 0, 'max': 1, 'palette': ['brown', 'yellow', 'green']},
             'NDVI')
    """)

    return {
        'standard': result_standard,
        'mowing': result_mowing,
        'difference': aoa_diff,
        'vis_params': {
            'di_standard': di_vis_standard,
            'di_mowing': di_vis_mowing,
            'aoa': aoa_vis,
            'diff': diff_vis
        }
    }


# =============================================================================
# STEP 7: Export Results
# =============================================================================

def export_results(result_mowing, result_standard, region):
    """
    Export AOA results to GEE assets for further analysis.
    """
    print("\n" + "="*70)
    print("EXPORTING RESULTS")
    print("="*70)

    # Export mowing-specific AOA (recommended)
    task_mowing = export_di_to_asset(
        di_image=result_mowing,
        asset_id='users/YOUR_USERNAME/aoa_mowing_alps_carpathians',  # UPDATE THIS
        region=region,
        scale=30,
        max_pixels=1e9
    )
    print(f"Mowing AOA export started: {task_mowing.id}")

    # Export standard AOA (for comparison)
    task_standard = export_di_to_asset(
        di_image=result_standard,
        asset_id='users/YOUR_USERNAME/aoa_standard_alps_carpathians',  # UPDATE THIS
        region=region,
        scale=30,
        max_pixels=1e9
    )
    print(f"Standard AOA export started: {task_standard.id}")

    print("\nMonitor exports at: https://code.earthengine.google.com/tasks")

    return task_mowing, task_standard


# =============================================================================
# STEP 8: Interpret Results for Mowing Classification
# =============================================================================

def interpret_results(train_di_mowing):
    """
    Provide interpretation guidelines for the mowing-specific AOA.
    """
    print("\n" + "="*70)
    print("INTERPRETATION GUIDE FOR MOWING-SPECIFIC AOA")
    print("="*70)

    print(f"""
Your mowing-specific AOA has been calculated with:
- Threshold: {train_di_mowing.threshold:.4f}
- Reference samples: {len(train_di_mowing.train_scaled)} mowing areas

HOW TO INTERPRET THE RESULTS:

1. DISSIMILARITY INDEX (DI):
   - Low DI (< {train_di_mowing.threshold:.4f}): Area is similar to known mowing areas
   - High DI (> {train_di_mowing.threshold:.4f}): Area is dissimilar to known mowing areas

2. AREA OF APPLICABILITY (AOA):
   - AOA = 1 (green): Model can reliably predict mowing in these areas
   - AOA = 0 (red): Predictions are uncertain - area is different from training mowing sites

3. FOR YOUR MOWING/NOMOWING CLASSIFICATION:

   a) Trust mowing predictions where:
      ✓ AOA = 1 (inside AOA)
      ✓ Model predicts "mowing"
      → High confidence mowing areas

   b) Be cautious where:
      ! AOA = 0 (outside AOA)
      ! Model predicts "mowing"
      → Low confidence - novel conditions, might be false positive

   c) Consider field validation where:
      - AOA = 0 but model predicts mowing with high probability
      - These are the most uncertain predictions

4. COMPARISON WITH SLOVAKIA/SWISS TEST DATA:
   - If your Slovakia/Swiss test areas show AOA = 1, it confirms they were
     within the training data's feature space
   - If AOA = 0 in some test areas but F1 scores were still good, this suggests
     your model generalizes well even outside the strict AOA

5. USING UNCERTAINTY FOR DECISION-MAKING:
   - Create "confidence zones":
     * High confidence: AOA = 1, DI < {train_di_mowing.threshold/2:.4f}
     * Medium confidence: AOA = 1, DI = {train_di_mowing.threshold/2:.4f}-{train_di_mowing.threshold:.4f}
     * Low confidence: AOA = 0

   - Use these zones to:
     * Prioritize field validation in low confidence areas
     * Weight mowing area statistics by confidence
     * Plan targeted data collection to improve model coverage

6. NEXT STEPS:
   - Compare DI values in Slovakia/Swiss test areas (should be inside AOA)
   - Check if high-elevation Alps have high DI (potential novel conditions)
   - Use DI to identify where additional training data would be most valuable
    """)


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    """
    Run complete AOA analysis for mowing classification.
    """
    print("="*70)
    print("AOA ANALYSIS FOR MOWING CLASSIFICATION")
    print("Alps and Carpathians Region")
    print("="*70)

    # Step 1: Load data
    print("\n1. Loading training data...")
    training_data = load_training_data()
    predictor_names = load_predictor_names()
    weights = load_variable_importance()
    cv_folds = create_cv_folds(training_data)

    print(f"  - Training samples: {len(training_data)}")
    print(f"  - Mowing samples: {sum(1 for s in training_data if s['class'] == 'mowing')}")
    print(f"  - NoMowing samples: {sum(1 for s in training_data if s['class'] == 'noMowing')}")
    print(f"  - Predictors: {', '.join(predictor_names)}")
    print(f"  - CV folds: {max(cv_folds) + 1}")

    # Create feature collection
    train_fc = create_training_fc(training_data, predictor_names)

    # Step 2: Calculate standard AOA (for comparison)
    train_di_standard = calculate_standard_aoa(train_fc, predictor_names, weights, cv_folds)

    # Step 3: Calculate mowing-specific AOA (recommended approach)
    train_di_mowing = calculate_mowing_specific_aoa(train_fc, predictor_names, weights, cv_folds)

    # Step 4: Load prediction area
    print("\n4. Loading prediction image (Alps and Carpathians)...")
    prediction_image, region = load_prediction_image()

    # Step 5: Apply both AOA approaches
    result_standard = apply_aoa_to_alps_carpathians(
        train_di_standard, prediction_image, region, "Standard AOA"
    )

    result_mowing = apply_aoa_to_alps_carpathians(
        train_di_mowing, prediction_image, region, "Mowing-specific AOA"
    )

    # Step 6: Compare and visualize
    comparison = compare_and_visualize(
        result_standard, result_mowing, region, prediction_image
    )

    # Step 7: Export (optional - uncomment to run)
    # export_results(result_mowing, result_standard, region)

    # Step 8: Interpretation guide
    interpret_results(train_di_mowing)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)

    return {
        'train_di_standard': train_di_standard,
        'train_di_mowing': train_di_mowing,
        'result_standard': result_standard,
        'result_mowing': result_mowing,
        'comparison': comparison,
        'prediction_image': prediction_image,
        'region': region
    }


# =============================================================================
# ADDITIONAL UTILITY: Check AOA for Test Data
# =============================================================================

def check_test_data_aoa(train_di_mowing, slovakia_swiss_data, predictor_names):
    """
    Check if your Slovakia/Swiss test data falls within the AOA.

    This helps validate that your good F1 scores on test data
    came from areas where the model should be reliable.
    """
    print("\n" + "="*70)
    print("CHECKING TEST DATA AOA")
    print("="*70)

    # Create test feature collection
    test_fc = create_training_fc(slovakia_swiss_data, predictor_names)

    # Calculate DI for test data
    test_list = test_fc.getInfo()['features']
    test_data = []
    for feat in test_list:
        props = feat['properties']
        row = [props.get(var, np.nan) for var in predictor_names]
        test_data.append(row)

    test_array = np.array(test_data, dtype=np.float64)

    # Scale test data
    test_scaled = np.zeros_like(test_array)
    for i, var in enumerate(predictor_names):
        test_scaled[:, i] = (test_array[:, i] - train_di_mowing.scale_mean[var]) / \
                            train_di_mowing.scale_std[var]
        test_scaled[:, i] *= train_di_mowing.weights[var]

    # Calculate DI for each test point
    test_dis = []
    for i in range(len(test_scaled)):
        # Find minimum distance to training data
        distances = np.sqrt(np.sum((train_di_mowing.train_scaled - test_scaled[i:i+1])**2, axis=1))
        min_dist = np.min(distances)
        di = min_dist / train_di_mowing.train_dist_avrgmean
        test_dis.append(di)

        aoa_status = "✓ Inside AOA" if di <= train_di_mowing.threshold else "✗ Outside AOA"
        print(f"  Test point {i+1}: DI = {di:.4f} {aoa_status}")

    # Summary statistics
    test_dis = np.array(test_dis)
    inside_aoa = np.sum(test_dis <= train_di_mowing.threshold)
    percentage = (inside_aoa / len(test_dis)) * 100

    print(f"\nTest Data Summary:")
    print(f"  - Total test points: {len(test_dis)}")
    print(f"  - Inside AOA: {inside_aoa} ({percentage:.1f}%)")
    print(f"  - Outside AOA: {len(test_dis) - inside_aoa} ({100-percentage:.1f}%)")
    print(f"  - Mean DI: {np.mean(test_dis):.4f}")
    print(f"  - Median DI: {np.median(test_dis):.4f}")

    return test_dis


# Run the analysis
if __name__ == "__main__":
    results = main()

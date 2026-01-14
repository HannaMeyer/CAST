"""
Area of Applicability (AOA) for Google Earth Engine
====================================================

This module implements the Area of Applicability (AOA) methodology from Meyer & Pebesma (2021)
for use with Google Earth Engine Python API.

Reference:
Meyer, H., Pebesma, E. (2021): Predicting into unknown space?
Estimating the area of applicability of spatial prediction models.
Methods in Ecology and Evolution 12: 1620-1633. https://doi.org/10.1111/2041-210X.13650

Author: Claude Code (based on CAST R package)
License: GPL (>= 2)
"""

import ee
import numpy as np
from typing import Dict, List, Union, Optional


class TrainDI:
    """
    Calculate Dissimilarity Index parameters from training data.

    This class computes the reference statistics needed to calculate AOA,
    including scaling parameters, distance statistics, and the DI threshold.

    Attributes:
        train_fc (ee.FeatureCollection): Training data as feature collection
        predictor_names (List[str]): List of predictor variable names
        weights (Dict[str, float]): Variable importance weights
        scale_mean (Dict[str, float]): Mean values for scaling
        scale_std (Dict[str, float]): Standard deviation values for scaling
        train_dist_avrgmean (float): Mean of average distances between training points
        threshold (float): DI threshold for determining AOA
        train_di (List[float]): DI values of training data
    """

    def __init__(self,
                 train: Union[ee.FeatureCollection, List[Dict]],
                 predictor_names: List[str],
                 weights: Optional[Dict[str, float]] = None,
                 cv_folds: Optional[List[int]] = None,
                 method: str = "L2",
                 use_cv: bool = False):
        """
        Initialize TrainDI and compute training statistics.

        Args:
            train: Training data as ee.FeatureCollection or list of dictionaries
            predictor_names: List of predictor variable names
            weights: Dictionary of variable importance weights (default: all 1.0)
            cv_folds: List of CV fold assignments for each training point (optional)
            method: Distance metric ("L2" for Euclidean, currently only L2 supported)
            use_cv: Whether to use CV folds for threshold calculation
        """
        # Convert to ee.FeatureCollection if needed
        if isinstance(train, list):
            features = [ee.Feature(None, feat) for feat in train]
            self.train_fc = ee.FeatureCollection(features)
        else:
            self.train_fc = train

        self.predictor_names = predictor_names
        self.method = method
        self.use_cv = use_cv
        self.cv_folds = cv_folds

        # Set weights (default to 1.0 for all variables)
        if weights is None:
            self.weights = {var: 1.0 for var in predictor_names}
        else:
            self.weights = weights

        # Compute training statistics
        self._compute_train_statistics()

    def _compute_train_statistics(self):
        """
        Compute scaling parameters and distance statistics from training data.

        This is the core algorithm that:
        1. Computes mean and std for z-score normalization
        2. Calculates average distances between all training points
        3. Computes minimum distances (considering CV folds if provided)
        4. Derives the DI threshold
        """
        # Step 1: Compute mean and standard deviation for each predictor
        # In GEE, we need to extract these from the FeatureCollection
        self.scale_mean = {}
        self.scale_std = {}

        # Get training data as a list (client-side operation)
        # Note: For large datasets, this should be done server-side in batches
        train_list = self.train_fc.getInfo()['features']

        # Convert to numpy array for efficient computation
        train_data = []
        for feat in train_list:
            props = feat['properties']
            row = [props.get(var, np.nan) for var in self.predictor_names]
            train_data.append(row)

        train_array = np.array(train_data, dtype=np.float64)
        n_samples = train_array.shape[0]

        # Compute scaling parameters
        for i, var in enumerate(self.predictor_names):
            col_data = train_array[:, i]
            self.scale_mean[var] = float(np.nanmean(col_data))
            self.scale_std[var] = float(np.nanstd(col_data, ddof=1))

        # Step 2: Scale and weight training data
        train_scaled = np.zeros_like(train_array)
        for i, var in enumerate(self.predictor_names):
            # Z-score normalization
            train_scaled[:, i] = (train_array[:, i] - self.scale_mean[var]) / self.scale_std[var]
            # Apply variable importance weights
            train_scaled[:, i] *= self.weights[var]

        # Step 3: Calculate distances between all training points
        # This implements the trainDI algorithm from CAST
        train_dist_avrg = []
        train_dist_min = []

        for i in range(n_samples):
            # Calculate distances from point i to all other points
            point_i = train_scaled[i:i+1, :]  # shape (1, n_predictors)

            # Euclidean distance to all points
            if self.method == "L2":
                distances = np.sqrt(np.sum((train_scaled - point_i)**2, axis=1))
            else:
                raise NotImplementedError(f"Method {self.method} not implemented")

            # Set distance to self as NA
            distances[i] = np.nan

            # Average distance (for normalization)
            train_dist_avrg.append(np.nanmean(distances))

            # Minimum distance (for DI calculation)
            # If using CV, mask out points not in the training fold
            if self.use_cv and self.cv_folds is not None:
                fold_i = self.cv_folds[i]
                # Mask distances to points in the same test fold
                dist_masked = distances.copy()
                for j in range(n_samples):
                    if self.cv_folds[j] == fold_i:
                        dist_masked[j] = np.nan

                if np.all(np.isnan(dist_masked)):
                    train_dist_min.append(np.nan)
                else:
                    train_dist_min.append(np.nanmin(dist_masked))
            else:
                train_dist_min.append(np.nanmin(distances))

        # Step 4: Calculate trainDist_avrgmean (normalization factor)
        self.train_dist_avrgmean = float(np.nanmean(train_dist_avrg))

        # Step 5: Calculate Training DI
        train_di = np.array(train_dist_min) / self.train_dist_avrgmean
        self.train_di = train_di[~np.isnan(train_di)].tolist()

        # Step 6: Calculate DI threshold (Q3 + 1.5 * IQR)
        q3 = float(np.nanpercentile(train_di, 75))
        iqr = float(np.nanpercentile(train_di, 75) - np.nanpercentile(train_di, 25))
        threshold = q3 + 1.5 * iqr

        # Cap threshold at maximum training DI
        max_train_di = float(np.nanmax(train_di))
        self.threshold = min(threshold, max_train_di)

        # Store scaled training data for distance calculations
        self.train_scaled = train_scaled

    def to_dict(self) -> Dict:
        """
        Export TrainDI parameters as a dictionary.

        Returns:
            Dictionary containing all parameters needed for AOA calculation
        """
        return {
            'predictor_names': self.predictor_names,
            'weights': self.weights,
            'scale_mean': self.scale_mean,
            'scale_std': self.scale_std,
            'train_dist_avrgmean': self.train_dist_avrgmean,
            'threshold': self.threshold,
            'method': self.method,
            'train_scaled': self.train_scaled.tolist()  # For distance calculations
        }

    @classmethod
    def from_dict(cls, params: Dict):
        """
        Create TrainDI object from parameter dictionary.

        Args:
            params: Dictionary containing TrainDI parameters

        Returns:
            TrainDI object with loaded parameters
        """
        # Create a dummy instance
        obj = cls.__new__(cls)
        obj.predictor_names = params['predictor_names']
        obj.weights = params['weights']
        obj.scale_mean = params['scale_mean']
        obj.scale_std = params['scale_std']
        obj.train_dist_avrgmean = params['train_dist_avrgmean']
        obj.threshold = params['threshold']
        obj.method = params['method']
        obj.train_scaled = np.array(params['train_scaled'])
        return obj


class AOA:
    """
    Area of Applicability calculator for GEE.

    This class applies the AOA methodology to prediction data in GEE,
    calculating the Dissimilarity Index (DI) and determining whether
    new locations are inside or outside the Area of Applicability.
    """

    def __init__(self, train_di: TrainDI):
        """
        Initialize AOA calculator with training statistics.

        Args:
            train_di: TrainDI object containing training statistics
        """
        self.train_di = train_di

    def calculate_di_image(self,
                          image: ee.Image,
                          scale: int = 30,
                          max_pixels: int = 1e8) -> ee.Image:
        """
        Calculate Dissimilarity Index for an ee.Image.

        This function applies the AOA algorithm to a multi-band image where
        each band represents a predictor variable.

        Args:
            image: ee.Image with bands matching predictor_names
            scale: Scale in meters for computation (default: 30m)
            max_pixels: Maximum number of pixels to process

        Returns:
            ee.Image with two bands: 'DI' and 'AOA'
        """
        # Select and order bands according to predictor names
        predictors = image.select(self.train_di.predictor_names)

        # Step 1: Scale predictors using training statistics
        scaled_bands = []
        for var in self.train_di.predictor_names:
            band = predictors.select([var])
            # Z-score normalization
            scaled = band.subtract(self.train_di.scale_mean[var]).divide(self.train_di.scale_std[var])
            # Apply variable importance weight
            weighted = scaled.multiply(self.train_di.weights[var])
            scaled_bands.append(weighted)

        # Combine scaled bands into single image
        scaled_image = ee.Image.cat(scaled_bands).rename(self.train_di.predictor_names)

        # Step 2: Calculate minimum distance to training data
        # Convert training data to ee.Array for efficient computation
        train_array = ee.Array(self.train_di.train_scaled.tolist())
        n_train = train_array.length().get([0])

        # Use a mapped approach to calculate distances
        def calc_min_distance(pixel_values):
            """Calculate minimum Euclidean distance from pixel to all training points."""
            # pixel_values is an ee.Array of shape [n_predictors]
            # Reshape to [1, n_predictors] for broadcasting
            pixel = ee.Array([pixel_values])

            # Calculate squared differences: (pixel - train_array)^2
            # This will broadcast to shape [n_train, n_predictors]
            diff = pixel.repeat(0, n_train).subtract(train_array)
            diff_squared = diff.pow(2)

            # Sum along predictor axis and take sqrt to get Euclidean distances
            distances = diff_squared.reduce(ee.Reducer.sum(), [1]).sqrt()

            # Return minimum distance
            return distances.reduce(ee.Reducer.min(), [0]).get([0])

        # Convert image to array image for efficient processing
        array_img = scaled_image.toArray()

        # Calculate minimum distance for each pixel
        min_dist_img = array_img.arrayReduce(ee.Reducer.toList(), [0]) \
            .arrayFlatten([['values']]) \
            .toArray() \
            .arrayProject([0]) \
            .arrayMap(calc_min_distance) \
            .arrayFlatten([['min_dist']])

        # Step 3: Calculate DI (normalized distance)
        di_image = min_dist_img.divide(self.train_di.train_dist_avrgmean).rename('DI')

        # Step 4: Calculate AOA (binary mask)
        aoa_image = di_image.lte(self.train_di.threshold).rename('AOA').toByte()

        # Combine DI and AOA into single image
        result = di_image.addBands(aoa_image)

        return result

    def calculate_di_optimized(self,
                              image: ee.Image,
                              geometry: ee.Geometry = None,
                              scale: int = 30,
                              tile_scale: int = 1) -> ee.Image:
        """
        Optimized DI calculation using neighborhoods and reduceRegion.

        This is a more GEE-native approach that handles large areas better.

        Args:
            image: ee.Image with predictor bands
            geometry: Region of interest (optional)
            scale: Scale in meters for computation
            tile_scale: Tile scale for large computations (1-16)

        Returns:
            ee.Image with 'DI' and 'AOA' bands
        """
        # Select and scale predictors
        scaled_bands = []
        for var in self.train_di.predictor_names:
            band = image.select([var])
            scaled = band.subtract(self.train_di.scale_mean[var]) \
                        .divide(self.train_di.scale_std[var]) \
                        .multiply(self.train_di.weights[var])
            scaled_bands.append(scaled)

        scaled_image = ee.Image.cat(scaled_bands)

        # Create training points image for distance calculation
        # This creates a multi-band image where each band represents one training point
        train_images = []
        for i in range(len(self.train_di.train_scaled)):
            train_point = self.train_di.train_scaled[i]

            # Create constant image for each training point
            point_bands = []
            for j, var in enumerate(self.train_di.predictor_names):
                point_bands.append(ee.Image.constant(train_point[j]))
            point_img = ee.Image.cat(point_bands)

            # Calculate Euclidean distance
            diff = scaled_image.subtract(point_img)
            dist = diff.pow(2).reduce(ee.Reducer.sum()).sqrt()
            train_images.append(dist)

        # Stack all distances and find minimum
        dist_stack = ee.Image.cat(train_images)
        min_dist = dist_stack.reduce(ee.Reducer.min())

        # Calculate DI
        di_image = min_dist.divide(self.train_di.train_dist_avrgmean).rename('DI')

        # Calculate AOA
        aoa_image = di_image.lte(self.train_di.threshold).rename('AOA').toByte()

        return di_image.addBands(aoa_image)


def calculate_aoa_for_class(train_fc: ee.FeatureCollection,
                            predictor_names: List[str],
                            class_column: str,
                            target_class: Union[str, int],
                            weights: Optional[Dict[str, float]] = None,
                            cv_folds: Optional[List[int]] = None) -> TrainDI:
    """
    Calculate TrainDI for a specific class only.

    This function addresses the class imbalance problem by computing
    AOA statistics based only on samples from the target class.

    Args:
        train_fc: Training feature collection
        predictor_names: List of predictor variable names
        class_column: Name of the column containing class labels
        target_class: The class value to use for AOA calculation (e.g., "mowing")
        weights: Variable importance weights (optional)
        cv_folds: CV fold assignments (optional)

    Returns:
        TrainDI object computed using only the target class samples

    Example:
        # Calculate AOA specifically for mowing areas
        >>> mowing_train_di = calculate_aoa_for_class(
        ...     train_fc=training_data,
        ...     predictor_names=['NDVI', 'elevation', 'slope'],
        ...     class_column='land_use',
        ...     target_class='mowing',
        ...     weights=variable_importance
        ... )
        >>>
        >>> # Apply to prediction area
        >>> aoa_calc = AOA(mowing_train_di)
        >>> result = aoa_calc.calculate_di_optimized(prediction_image)
    """
    # Filter training data to target class only
    class_filtered = train_fc.filter(ee.Filter.eq(class_column, target_class))

    # Get filtered CV folds if provided
    filtered_cv_folds = None
    if cv_folds is not None:
        # Get indices of samples in target class
        train_list = train_fc.getInfo()['features']
        target_indices = [i for i, feat in enumerate(train_list)
                         if feat['properties'].get(class_column) == target_class]
        filtered_cv_folds = [cv_folds[i] for i in target_indices]

    # Calculate TrainDI using only target class
    train_di = TrainDI(
        train=class_filtered,
        predictor_names=predictor_names,
        weights=weights,
        cv_folds=filtered_cv_folds,
        use_cv=(filtered_cv_folds is not None)
    )

    return train_di


def create_training_fc(data: List[Dict],
                      predictor_names: List[str],
                      lat_col: str = 'latitude',
                      lon_col: str = 'longitude') -> ee.FeatureCollection:
    """
    Create an ee.FeatureCollection from training data.

    Args:
        data: List of dictionaries containing training samples
        predictor_names: List of predictor variable names
        lat_col: Column name for latitude (default: 'latitude')
        lon_col: Column name for longitude (default: 'longitude')

    Returns:
        ee.FeatureCollection with training data

    Example:
        >>> training_data = [
        ...     {'longitude': 10.5, 'latitude': 47.5, 'NDVI': 0.8, 'elevation': 1500, 'class': 'mowing'},
        ...     {'longitude': 10.6, 'latitude': 47.6, 'NDVI': 0.3, 'elevation': 800, 'class': 'noMowing'},
        ... ]
        >>> train_fc = create_training_fc(training_data, ['NDVI', 'elevation'])
    """
    features = []
    for sample in data:
        # Extract coordinates
        lon = sample.get(lon_col)
        lat = sample.get(lat_col)

        # Create geometry (if coordinates provided)
        if lon is not None and lat is not None:
            geom = ee.Geometry.Point([lon, lat])
        else:
            geom = None

        # Create properties dict with predictors and other attributes
        properties = {k: v for k, v in sample.items()}

        # Create feature
        features.append(ee.Feature(geom, properties))

    return ee.FeatureCollection(features)


# Utility functions for working with GEE AOA

def export_di_to_asset(di_image: ee.Image,
                       asset_id: str,
                       region: ee.Geometry,
                       scale: int = 30,
                       max_pixels: int = 1e9):
    """
    Export DI/AOA image to GEE asset.

    Args:
        di_image: Image with DI and AOA bands
        asset_id: Asset ID for export (e.g., 'users/username/aoa_alps')
        region: Region to export
        scale: Export scale in meters
        max_pixels: Maximum number of pixels to export
    """
    task = ee.batch.Export.image.toAsset(
        image=di_image,
        description='AOA_export',
        assetId=asset_id,
        region=region,
        scale=scale,
        maxPixels=max_pixels,
        crs='EPSG:4326'
    )
    task.start()
    return task


def visualize_aoa(di_image: ee.Image,
                 region: ee.Geometry = None,
                 scale: int = 1000) -> Dict:
    """
    Create visualization parameters for AOA results.

    Args:
        di_image: Image with DI and AOA bands
        region: Region for statistics (optional)
        scale: Scale for statistics computation

    Returns:
        Dictionary with visualization parameters
    """
    # Calculate DI statistics
    if region is not None:
        stats = di_image.select('DI').reduceRegion(
            reducer=ee.Reducer.percentile([5, 50, 95]),
            geometry=region,
            scale=scale,
            maxPixels=1e9
        ).getInfo()

        p5 = stats.get('DI_p5', 0)
        p95 = stats.get('DI_p95', 1)
    else:
        p5, p95 = 0, 1

    vis_params = {
        'DI': {
            'bands': ['DI'],
            'min': p5,
            'max': p95,
            'palette': ['blue', 'white', 'red']
        },
        'AOA': {
            'bands': ['AOA'],
            'min': 0,
            'max': 1,
            'palette': ['red', 'green']
        }
    }

    return vis_params


# Example usage and documentation
__doc__ += """

Example Usage:
==============

1. Basic AOA calculation:
--------------------------

```python
import ee
from gee_aoa import TrainDI, AOA, create_training_fc

# Initialize Earth Engine
ee.Initialize()

# Prepare training data
training_data = [
    {'longitude': 10.5, 'latitude': 47.5, 'NDVI': 0.8, 'elevation': 1500},
    {'longitude': 10.6, 'latitude': 47.6, 'NDVI': 0.3, 'elevation': 800},
    # ... more samples
]

train_fc = create_training_fc(training_data, ['NDVI', 'elevation'])

# Calculate training DI
train_di = TrainDI(
    train=train_fc,
    predictor_names=['NDVI', 'elevation'],
    weights={'NDVI': 0.7, 'elevation': 0.3}
)

# Load prediction image
prediction_image = ee.Image('your_composite_image')

# Calculate AOA
aoa_calculator = AOA(train_di)
result = aoa_calculator.calculate_di_optimized(prediction_image, scale=30)

# Visualize
print('DI Threshold:', train_di.threshold)
Map.addLayer(result.select('DI'), {'min': 0, 'max': 2, 'palette': ['blue', 'white', 'red']}, 'DI')
Map.addLayer(result.select('AOA'), {'min': 0, 'max': 1, 'palette': ['red', 'green']}, 'AOA')
```

2. Class-specific AOA for mowing/noMowing:
-------------------------------------------

```python
from gee_aoa import calculate_aoa_for_class

# Training data with class labels
training_data_with_class = [
    {'longitude': 10.5, 'latitude': 47.5, 'NDVI': 0.8, 'elevation': 1500, 'class': 'mowing'},
    {'longitude': 10.6, 'latitude': 47.6, 'NDVI': 0.3, 'elevation': 800, 'class': 'noMowing'},
    # ... many more noMowing samples
]

train_fc = create_training_fc(training_data_with_class, ['NDVI', 'elevation'])

# Calculate AOA ONLY for mowing class
mowing_train_di = calculate_aoa_for_class(
    train_fc=train_fc,
    predictor_names=['NDVI', 'elevation'],
    class_column='class',
    target_class='mowing',  # Only use mowing samples as reference
    weights={'NDVI': 0.7, 'elevation': 0.3}
)

# Apply to prediction area
aoa_calc = AOA(mowing_train_di)
result = aoa_calc.calculate_di_optimized(alps_carpathians_image, scale=30)

# Now DI shows dissimilarity to MOWING areas specifically
# AOA shows where predictions are reliable for mowing classification
Map.addLayer(result.select('AOA'), {'min': 0, 'max': 1, 'palette': ['red', 'green']},
             'AOA for Mowing')
```

3. Using CV folds:
------------------

```python
# With spatial cross-validation folds
cv_folds = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]  # Region IDs

train_di = TrainDI(
    train=train_fc,
    predictor_names=['NDVI', 'elevation', 'slope'],
    weights=variable_importance,
    cv_folds=cv_folds,
    use_cv=True  # Use CV folds for threshold calculation
)
```

4. Export results:
------------------

```python
from gee_aoa import export_di_to_asset

# Define region (Alps and Carpathians)
region = ee.Geometry.Rectangle([5, 43, 27, 49])

# Export to asset
task = export_di_to_asset(
    di_image=result,
    asset_id='users/your_username/aoa_mowing_alps',
    region=region,
    scale=30
)

# Monitor task
print(f'Task started: {task.id}')
```

Key Advantages for Mowing Classification:
==========================================

1. **Class-specific AOA**: By using only mowing samples, the AOA shows where
   your model can reliably predict mowing areas, not just where it can predict
   in general.

2. **Handles class imbalance**: With many noMowing samples, standard AOA would
   be dominated by noMowing. Class-specific AOA focuses on the minority class.

3. **Interpretability**: Areas inside the mowing-specific AOA are similar to
   known mowing areas in the predictor space, giving you confidence in mowing
   predictions there.

4. **Uncertainty quantification**: DI values give you continuous uncertainty
   estimates - higher DI means less similar to training mowing areas.
"""

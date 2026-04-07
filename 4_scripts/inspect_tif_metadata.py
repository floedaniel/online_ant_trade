#!/usr/bin/env python3
"""Inspect GeoTIFF metadata for current and future bioclim files"""

import os
try:
    from osgeo import gdal
    gdal.UseExceptions()

    print("\n=== INSPECTING CURRENT AND FUTURE BIOCLIM RASTERS ===\n")

    # File paths
    current_file = "./1_raw_data/bioclim/current/current.tif"
    future_file = "./1_raw_data/bioclim/future/future.tif"

    for label, filepath in [("CURRENT", current_file), ("FUTURE", future_file)]:
        print(f"\n{'='*70}")
        print(f"=== {label} BIOCLIM ===")
        print(f"File: {filepath}")

        if not os.path.exists(filepath):
            print(f"ERROR: File not found: {filepath}")
            continue

        # Open dataset
        ds = gdal.Open(filepath)
        if ds is None:
            print(f"ERROR: Could not open {filepath}")
            continue

        # Basic info
        print(f"\nNumber of bands/layers: {ds.RasterCount}")
        print(f"Dimensions (width x height): {ds.RasterXSize} x {ds.RasterYSize}")
        print(f"Total cells per layer: {ds.RasterXSize * ds.RasterYSize:,}")

        # Geotransform
        gt = ds.GetGeoTransform()
        print(f"\nGeotransform:")
        print(f"  Origin (upper left): ({gt[0]:.6f}, {gt[3]:.6f})")
        print(f"  Pixel size: {gt[1]:.6f} x {gt[5]:.6f}")

        # Extent
        xmin = gt[0]
        xmax = gt[0] + ds.RasterXSize * gt[1]
        ymax = gt[3]
        ymin = gt[3] + ds.RasterYSize * gt[5]
        print(f"\nExtent:")
        print(f"  X: {xmin:.6f} to {xmax:.6f}")
        print(f"  Y: {ymin:.6f} to {ymax:.6f}")

        # CRS
        proj = ds.GetProjection()
        if proj:
            print(f"\nProjection: {proj[:100]}...")

        # Band/layer information
        print(f"\nBand/Layer names and statistics:")
        for i in range(1, min(ds.RasterCount + 1, 20)):  # Show first 20
            band = ds.GetRasterBand(i)
            band_name = band.GetDescription()
            if not band_name:
                band_name = f"Band_{i}"

            # Get statistics if available
            stats = band.GetStatistics(False, False)
            if stats:
                print(f"  {i:2d}. {band_name:40s} | Min: {stats[0]:10.2f} | Max: {stats[1]:10.2f} | Mean: {stats[2]:10.2f}")
            else:
                print(f"  {i:2d}. {band_name}")

        if ds.RasterCount > 20:
            print(f"  ... and {ds.RasterCount - 20} more layers")

        ds = None  # Close

    print(f"\n{'='*70}")
    print("\n=== INSPECTION COMPLETE ===\n")

except ImportError:
    print("GDAL Python bindings not available. Trying alternative method...")

    # Alternative: just show basic file info
    import os
    current_file = "./1_raw_data/bioclim/current/current.tif"
    future_file = "./1_raw_data/bioclim/future/future.tif"

    print("\n=== BASIC FILE INFORMATION ===\n")
    for label, filepath in [("CURRENT", current_file), ("FUTURE", future_file)]:
        if os.path.exists(filepath):
            size_mb = os.path.getsize(filepath) / (1024 * 1024)
            print(f"{label}: {filepath}")
            print(f"  File size: {size_mb:.1f} MB")
            print()

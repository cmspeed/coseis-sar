import xarray as xr
import numpy as np
from pathlib import Path

wavelength = 0.056

def convert_file(in_path, out_path):
    ds = xr.open_dataset(in_path)

    # Create new variable
    los_disp = ds["losDisplacement"]
    unwrapped_phase = (los_disp * (4 * np.pi / wavelength)).astype("float32")

    # Add it as a new variable
    ds["unwrappedPhase"] = unwrapped_phase
    ds["unwrappedPhase"].attrs["units"] = "radians"

    # Remove losDisplacement
    ds = ds.drop_vars("losDisplacement")

    # Save new file (overwrite or new path)
    ds.to_netcdf(out_path)

# Example: process all .nc files in a directory
input_dir = Path("data/")
for nc in input_dir.glob("*.nc"):
    convert_file(nc, nc)  # overwrite in place

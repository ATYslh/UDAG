import json
import os
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import regionmask
import xarray as xr
from cdo import Cdo


def mask_data(input_file):
    world = gpd.read_file("ne_10m_admin_0_countries.shp")
    germany = world[world["NAME"] == "Germany"].to_crs("EPSG:4326")

    ds = xr.open_dataset(input_file)

    # Create regionmask
    mask = regionmask.Regions(
        [germany.geometry.values[0]], names=["Germany"], abbrevs=["DE"]
    )
    germany_mask = mask.mask(ds.lon, ds.lat)

    # Drop vertices variables (not needed)
    vars_to_drop = [v for v in ds.data_vars if "vertices" in v]
    ds = ds.drop_vars(vars_to_drop)

    # Mask only variables that depend on time, rlat, and rlon
    spatial_vars = [
        v for v in ds.data_vars if {"time", "rlat", "rlon"} <= set(ds[v].dims)
    ]
    # Expand mask to have a time dimension
    mask3d = germany_mask.expand_dims(time=ds.time)

    # Apply only to variables with time, rlat, rlon
    masked = ds[["sfcWind"]].where(mask3d == 0, drop=True)

    # Merge back with untouched non-spatial variables
    ds_germany = xr.merge([masked, ds.drop_vars(spatial_vars)])

    # Drop boundary variables except time_bnds
    vars_to_drop = [
        v for v in ds_germany.data_vars if "bnds" in v and v not in ["time_bnds"]
    ]
    ds_germany = ds_germany.drop_vars(vars_to_drop)

    ds_germany.to_netcdf("/scratch/g/g260190/germany_masked_large.nc")
    return "/scratch/g/g260190/germany_masked_large.nc"


def write_json_file(filename: str, content: dict) -> None:
    with open(filename, "w") as file:
        json.dump(content, file, indent=4)


def find_folders(base, temporal_resolution) -> list[str]:
    if temporal_resolution == "yearly":
        temporal_resolution = "mon"
    matches = base.glob(
        f"CLMcom-*/*/*/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/{temporal_resolution}/sfcWind/v20240920"
    )
    list_folders = []
    for m in matches:
        if m.is_dir():
            list_folders.append(str(m))
    return sorted(list_folders)


def get_sorted_nc_files(folder_path: str, substring=None):
    """
    Returns a sorted list of all nc_files in the folder_path.
    """
    nc_files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(".nc")
        and os.path.isfile(os.path.join(folder_path, f))
        and (substring is None or substring in f)
    ]
    return sorted(nc_files)


def generate_filename(folder: str, variable: str) -> str:
    parts = folder.split("/")
    return "_".join([parts[6], parts[7], parts[8], parts[9], variable]) + ".nc"


def create_yearly_data(
    input_folder, output_folder, overwrite, temporal_resolution
) -> None:
    output_folder = os.path.join(output_folder, temporal_resolution)
    output_filename = os.path.join(
        output_folder, generate_filename(input_folder, "sfcWind")
    )
    if os.path.exists(output_filename) and not overwrite:
        return

    files = get_sorted_nc_files(input_folder)
    dummy_data = "/scratch/g/g260190/dummy.nc"
    print(f"Working on {input_folder} {datetime.now()}", file=sys.stderr)

    cdo = Cdo()

    # Apply mask and fldmean to each file individually before merging
    fldmean_files = []
    for f in files:
        # First mask the file
        masked_file = os.path.join(
            "/scratch/g/g260190/", os.path.basename(f).replace(".nc", "_masked.nc")
        )
        masked_file = mask_data(f)  # assuming mask_data returns a filename

        # Then compute fldmean on the masked file
        fldmean_file = os.path.join(
            "/scratch/g/g260190/", os.path.basename(f).replace(".nc", "_fldmean.nc")
        )
        cdo.fldmean(input=masked_file, output=fldmean_file)

        fldmean_files.append(fldmean_file)

    # Now merge the reduced files in time
    cdo.mergetime(input=fldmean_files, output=dummy_data)

    # Check if we have to limit the years to 2100
    ds = xr.open_dataset(dummy_data)
    years = ds["time"].dt.year.values
    max_year = years.max()
    if max_year >= 2100:
        limited_data = "dummy_data_2015_2100.nc"
        cdo.selyear("2015/2100", input=dummy_data, output=limited_data)
        os.system(f"mv {limited_data} {dummy_data}")

    # At this point dummy_data is already masked+fldmean, so just do temporal averaging
    if temporal_resolution == "yearly":
        cdo.yearmonmean(input=dummy_data, output=output_filename)
    else:
        os.system(f"mv {dummy_data} {output_filename}")


def create_info_json(output_folder):
    info = {}
    for file in get_sorted_nc_files(output_folder):
        parts = os.path.basename(file).split("_")
        temp_resolution = os.path.dirname(file).split("/")[-1]
        info.setdefault(temp_resolution, {}).setdefault(parts[2], {}).setdefault(
            parts[3], {}
        )[parts[0]] = file

    write_json_file("sfc_Wind_info.json", info)


if __name__ == "__main__":
    output_folder = "/work/bb1203/g260190_heinrich/UDAG/Data/sfcWind"
    overwrite = False

    for spatial_resolution in ["EUR-12", "MEU-3"]:
        for temporal_resolution in ["yearly", "mon", "day", "1hr"]:
            base = Path(f"/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/{spatial_resolution}")
            sfc_folders = find_folders(base, temporal_resolution)
            for input_folder in sfc_folders:
                create_yearly_data(
                    input_folder, output_folder, overwrite, temporal_resolution
                )

    create_info_json(output_folder)

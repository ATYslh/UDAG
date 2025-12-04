import json
import os
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import regionmask
import xarray as xr
from cdo import Cdo


def mask_data(input_file: str, mask2d, variable: str):
    """
    Compute Germany fldmean for variable and write a compact time series file.
    The output variable name remains 'variable'.
    """
    output_file = os.path.join(
        "/scratch/g/g260190/",
        os.path.basename(input_file).replace(".nc", "_fldmean.nc"),
    )

    # Open only the needed variable with dask chunks
    ds = xr.open_dataset(input_file, chunks={"time": 12})[[variable]]

    # Broadcast mask2d across time (no expand_dims needed)
    fldmean = ds[variable].where(mask2d == 0).mean(dim=["rlat", "rlon"], skipna=True)

    # Keep the variable name 'variable'
    ds_out = fldmean.to_dataset(name=variable)

    # Efficient write
    encoding = {variable: {"zlib": True, "complevel": 4}}
    ds_out.to_netcdf(output_file, engine="h5netcdf", encoding=encoding)

    return output_file


def write_json_file(filename: str, content: dict) -> None:
    with open(filename, "w") as file:
        json.dump(content, file, indent=4)


def find_folders(
    base, temporal_resolution, variable, project, spatial_resolution
) -> list[str]:
    if temporal_resolution == "yearly":
        temporal_resolution = "mon"

    if project == "UDAG":
        base = Path(f"/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/{spatial_resolution}")
        matches = base.glob(
            f"CLMcom-*/*/*/*/ICON-CLM-202407-1-1/v1-r1/{temporal_resolution}/{variable}/v20240920"
        )
    elif project == "NUKLEUS":
        base = Path("/work/bb1203/data_NUKLEUS_CMOR/CEU-3/")
        matches = base.glob(f"CLMcom-*/*/*/*/*/*/{temporal_resolution}/{variable}/*")
    else:
        raise ValueError(f"Unknown project: {project}")

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
    input_folder, output_folder, overwrite, temporal_resolution, variable: str
) -> None:
    output_folder = os.path.join(output_folder, temporal_resolution)
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_filename = os.path.join(
        output_folder, generate_filename(input_folder, variable)
    )
    if os.path.exists(output_filename) and not overwrite:
        return

    files = get_sorted_nc_files(input_folder)
    dummy_data = "/scratch/g/g260190/dummy.nc"
    print(f"Working on {input_folder} {datetime.now()}", file=sys.stderr)

    cdo = Cdo()

    # Apply mask and fldmean to each file individually before merging
    fldmean_files = []
    if "EUR-12" in input_folder:
        mask_file = "mask/mask_EUR-12.nc"
    elif "MEU-3" in input_folder:
        mask_file = "mask/mask_MEU-3.nc"
    elif "CEU-3" in input_folder:
        mask_file = "mas/mask_CEU-3.nc"
    else:
        raise ValueError(f"Unknown resolution in filename: {input_folder}")

    with xr.open_dataarray(mask_file) as mask2d:
        for f in files:
            fldmean_files.append(mask_data(f, mask2d, variable))

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


def sort_dict_recursively(d):
    """
    Recursively sorts a dictionary by its keys alphabetically.
    """
    if not isinstance(d, dict):
        return d  # Base case: return non-dict values as-is

    # Sort keys and apply recursively to values
    return {k: sort_dict_recursively(d[k]) for k in sorted(d)}


def create_info_json(output_folder):
    info = {}
    variable = output_folder.split("/")[-1]
    country = output_folder.split("/")[-2]
    project = output_folder.split("/")[-3]
    folders = [
        f
        for f in os.listdir(output_folder)
        if os.path.isdir(os.path.join(output_folder, f))
    ]
    for folder in folders:
        for file in get_sorted_nc_files(folder):
            parts = os.path.basename(file).split("_")
            temp_resolution = os.path.dirname(file).split("/")[-1]
            info.setdefault(temp_resolution, {}).setdefault(parts[2], {}).setdefault(
                parts[3], {}
            )[parts[0]] = file

    write_json_file(
        f"{project}_{country}_{variable}_info.json", sort_dict_recursively(info)
    )


def precompute_masks(country):
    resolutions = {
        "EUR-12": "/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/EUR-12/CLMcom-Hereon/EC-Earth3-Veg/historical/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/mon/sfcWind/v20240920/sfcWind_EUR-12_EC-Earth3-Veg_historical_r1i1p1f1_CLMcom-Hereon_ICON-CLM-202407-1-1_v1-r1_mon_195001-195012.nc",
        "MEU-3": "/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/MEU-3/CLMcom-Hereon/EC-Earth3-Veg/historical/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/mon/sfcWind/v20240920/sfcWind_MEU-3_EC-Earth3-Veg_historical_r1i1p1f1_CLMcom-Hereon_ICON-CLM-202407-1-1_v1-r1_mon_196001-196012.nc",
        "CEU-3": "/work/bb1203/data_NUKLEUS_CMOR/CEU-3/CLMcom-Hereon/MPI-M-MPI-ESM1-2-HR/historical/r1i1p1f1/CLMcom-Hereon-CCLM-6-0-clm2/nukleus-x2yn2-v1/mon/tas/v20230222/tas_CEU-3_MPI-M-MPI-ESM1-2-HR_historical_r1i1p1f1_CLMcom-Hereon-CCLM-6-0-clm2_nukleus-x2yn2-v1_mon_196101-196112.nc",
    }
    shapefile = "/work/bb1203/g260190_heinrich/UDAG/Scripts/shape_files/ne_10m_admin_0_countries.shp"

    world = gpd.read_file(shapefile)
    if country == "Germany":
        abbrev = "DE"
    elif country == "Denmark":
        abbrev = "DK"
    else:
        raise ValueError(f"Unknwon country {country}")

    country_info = world.loc[world["NAME"] == country].to_crs("EPSG:4326")
    region = regionmask.Regions(
        [country_info.geometry.values[0]], names=[country], abbrevs=[abbrev]
    )

    saved_masks = {}

    for res, nc_file in resolutions.items():
        ds = xr.open_dataset(nc_file)

        # Build mask for this resolution grid (2D only)
        mask2d = region.mask(ds.lon, ds.lat)

        # Save mask
        mask_file = f"masks/mask_{res}.nc"
        mask2d.to_netcdf(mask_file)

        saved_masks[res] = mask_file
        print(f"Saved mask for {res} to {mask_file}")

    return saved_masks


def main():
    variable = "pr"
    country = "Germany"

    project = "NUKLEUS"
    output_folder = (
        f"/work/bb1203/g260190_heinrich/UDAG/Data/{project}/{country}/{variable}"
    )
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    list_of_wanted_resolutions = ["yearly"]  # ["yearly", "mon", "day", "1hr"]

    overwrite = False
    precompute_masks(country)

    for spatial_resolution in ["EUR-12", "MEU-3", "CEU-3"]:
        if (spatial_resolution == "CEU-3" and project != "NUKLEUS") or (
            spatial_resolution != "CEU-3" and project == "NUKLEUS"
        ):
            continue
        for temporal_resolution in list_of_wanted_resolutions:
            data_folders = find_folders(
                project, temporal_resolution, variable, project, spatial_resolution
            )
            for input_folder in data_folders:
                create_yearly_data(
                    input_folder,
                    output_folder,
                    overwrite,
                    temporal_resolution,
                    variable,
                )

    create_info_json(output_folder)


if __name__ == "__main__":
    main()

import hashlib
import json
import os
import sys
from datetime import datetime

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
    if variable=="sfcWind":
        variable="fg"
    # Open only the needed variable with dask chunks
    ds = xr.open_dataset(input_file, chunks={})[[variable]]

    # Broadcast mask2d across time (no expand_dims needed)
    fldmean = ds[variable].where(mask2d == 0).mean(dim=["lat", "lon"], skipna=True)

    # Keep the variable name 'variable'
    if variable=="fg":
        variable="sfcWind"
    ds_out = fldmean.to_dataset(name=variable)

    encoding = {variable: {"zlib": True, "complevel": 4}}
    ds_out.to_netcdf(output_file, engine="h5netcdf", encoding=encoding)
    os.system(f"cdo settaxis,1980-01-01,00:00:00,1day {output_file} {output_file}z")
    os.system(f"mv {output_file}z {output_file}")
    return output_file


def write_json_file(filename: str, content: dict) -> None:
    with open(filename, "w") as file:
        json.dump(content, file, indent=4)


def get_sorted_finished_nc_files(folder_path: str, substring=None):
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


def get_sorted_nc_files(folder_path: str, substring=None):
    """
    Returns a sorted list of all nc_files in the folder_path.
    """

    return [
        "/work/gg0302/observations/eObs/eObs31.0/daily/0.1deg_reg/fg_ens_mean_0.1deg_reg_v31.0e.nc"
    ]


def create_yearly_data(
    output_folder, overwrite, list_temporal_resolution, variable: str
) -> None:

    files = [
        "/work/gg0302/observations/eObs/eObs31.0/daily/0.1deg_reg/fg_ens_mean_0.1deg_reg_v31.0e.nc"
    ]
    now = datetime.now().isoformat()
    hash_value = hashlib.sha256(now.encode()).hexdigest()

    print(
        f"Working on /work/gg0302/observations/eObs/eObs31.0/daily/0.1deg_reg/fg_ens_mean_0.1deg_reg_v31.0e.nc {datetime.now()}",
        file=sys.stderr,
    )

    cdo = Cdo()
    temporal_resolution=list_temporal_resolution[0]
    # Apply mask and fldmean to each file individually before merging
    fldmean_files = []
    if "eObs" in files[0]:
        mask_file = "/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/mask_EOBS.nc"
    else:
        raise ValueError(f"Unknown resolution in filename: {files[0]}")

    with xr.open_dataarray(mask_file) as mask2d:
        for f in files:
            fldmean_files.append(mask_data(f, mask2d, variable))
    dummy_data = fldmean_files[0]

    # Check if we have to limit the years to 2100
    ds = xr.open_dataset(dummy_data)
    years = ds["time"].dt.year.values
    max_year = years.max()
    if max_year >= 2100:
        limited_data = "dummy_data_2015_2100.nc"
        cdo.selyear("2015/2100", input=dummy_data, output=limited_data)
        os.system(f"mv {limited_data} {dummy_data}")


    # At this point dummy_data is already masked+fldmean, so just do temporal averaging
    print(list_temporal_resolution)
    for temporal_resolution in list_temporal_resolution:
        output_folder2 = os.path.join(output_folder, temporal_resolution)

        if not os.path.exists(output_folder2):
            os.mkdir(output_folder2)

        output_filename = os.path.join(output_folder2, "EOBS.nc")
        if os.path.exists(output_filename) and not overwrite:
            return
        print(f"output folder: {output_folder2} for {temporal_resolution}")
        if temporal_resolution == "yearly":
            cdo.yearmean(input=dummy_data, output=output_filename)
        elif temporal_resolution == "mon":
            cdo.monmean(input=dummy_data, output=output_filename)
        elif temporal_resolution == "day":
            os.system(f"cp {dummy_data} {output_filename}")
        else:
            raise ValueError("Unknown resolution or resolution not available.")

def sorted_resolution(list_of_wanted_resolutions) -> list[str]:
    temporary_resolutions = {
        "yearly": 1,
        "mon": 2,
        "day": 3,
        "1hr": 4,
    }
    
    return sorted(list_of_wanted_resolutions, key=lambda x: temporary_resolutions[x], reverse=True)

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
        os.path.join(output_folder, f)
        for f in os.listdir(output_folder)
        if os.path.isdir(os.path.join(output_folder, f))
    ]
    for folder in folders:
        for file in get_sorted_finished_nc_files(folder):
            temp_resolution = os.path.dirname(file).split("/")[-1]
            info[temp_resolution] = file

    write_json_file(
        f"/work/bb1364/g260190_heinrich/UDAG/Data/json_files/{project}_{country}_{variable}_info.json",
        sort_dict_recursively(info),
    )


def precompute_masks(country):
    resolutions = {
        "EOBS": "/work/gg0302/observations/eObs/eObs31.0/daily/0.1deg_reg/fg_ens_mean_0.1deg_reg_v31.0e.nc",
    }
    shapefile = "/work/bb1364/g260190_heinrich/UDAG/Scripts/shape_files/ne_10m_admin_0_countries.shp"

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
        mask_file = f"/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/mask_{res}.nc"
        mask2d.to_netcdf(mask_file)

        saved_masks[res] = mask_file
        print(f"Saved mask for {res} to {mask_file}")

    return saved_masks


def main():
    variables = ["sfcWind"]
    country = "Germany"
    project = "EOBS"

    list_of_wanted_resolutions = [
        "yearly",
        "mon",
        "day",
    ]  # ["yearly", "mon", "day", "1hr"]
    list_of_wanted_resolutions=sorted_resolution(list_of_wanted_resolutions)
    overwrite = True
    precompute_masks(country)

    for variable in variables:
        output_folder = (
            f"/work/bb1364/g260190_heinrich/UDAG/Data/{project}/{country}/{variable}"
        )

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)


        create_yearly_data(
            output_folder,
            overwrite,
            list_of_wanted_resolutions,
            variable,
        )

        create_info_json(output_folder)


if __name__ == "__main__":
    main()

import hashlib
import json
import os
import sys
from datetime import datetime

import geopandas as gpd
import regionmask
import xarray as xr
from cdo import Cdo

# TODO: Add rsds
# TODO: RENAME: sfcWind
# TODO REPLACE MV with copy for daily
# TODO: FIX Bullshit for frequency


def mask_data(input_file: str, mask2d, variable: str):
    """
    Compute Germany fldmean for variable and write a compact time series file.
    The output variable name remains 'variable'.
    """
    if variable == "tas":
        variable_old = "T_2M"
    elif variable == "tasmax":
        variable_old = "TMAX_2M"
    elif variable == "tasmin":
        variable_old = "TMIN_2M"
    elif variable == "pr":
        variable_old = "TOT_PREC"
    elif variable == "rsds":
        variable_old = "ASOD_S"
    elif variable == "sfcWind":
        variable_old = "SP_10M"
    else:
        raise ValueError("unknown variable")
    output_file = os.path.join(
        "/scratch/g/g260190/",
        os.path.basename(input_file).replace(
            ".nc", f"_{hashlib.md5(input_file.encode()).hexdigest()}.nc"
        ),
    )

    # Open only the needed variable with dask chunks
    ds = xr.open_dataset(input_file, chunks={})[[variable_old]]

    # Broadcast mask2d across time (no expand_dims needed)
    fldmean = (
        ds[variable_old].where(mask2d == 0).mean(dim=["rlat", "rlon"], skipna=True)
    )

    # Keep the variable name 'variable'
    ds_out = fldmean.to_dataset(name=variable_old)

    encoding = {variable_old: {"zlib": True, "complevel": 4}}
    ds_out.to_netcdf(output_file, engine="h5netcdf", encoding=encoding)

    return output_file


def write_json_file(filename: str, content: dict) -> None:
    with open(filename, "w") as file:
        json.dump(content, file, indent=4)


def find_folders(variable) -> list[str]:
    base = "/work/bb1364/production_runs/work/IAEVALH01/post/yearly"
    variable_old = ""

    if variable == "tas":
        variable_old = "T_2M"
    elif variable == "tasmax":
        variable_old = "TMAX_2M"
    elif variable == "tasmin":
        variable_old = "TMIN_2M"
    elif variable == "pr":
        variable_old = "TOT_PREC"
    elif variable == "rsds":
        variable_old = "ASOD_S"
    elif variable == "sfcWind":
        variable_old = "SP_10M"
    else:
        raise ValueError("unknown variable")
    return [os.path.join(base, variable_old)]


def get_sorted_nc_files(folder_path: str, substring=None):
    """
    Returns a sorted list of all nc_files in the folder_path.
    """
    nc_files = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.endswith(".ncz")
        and os.path.isfile(os.path.join(folder_path, f))
        and (substring is None or substring in f)
    ]
    return sorted(nc_files)


def generate_filename(variable: str) -> str:
    return f"MEU-3_CLMcom-Hereon_ERA5_evaluation_{variable}.nc"


def create_yearly_data(
    input_folder, output_folder, overwrite, temporal_resolution, variable: str
) -> None:
    output_folder = os.path.join(output_folder, temporal_resolution[0])
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_filename = os.path.join(output_folder, generate_filename(variable))
    if os.path.exists(output_filename) and not overwrite:
        return

    files = get_sorted_nc_files(input_folder)
    dummy_data = "/scratch/g/g260190/dummy.nc"
    print(f"Working on {input_folder} {datetime.now()}", file=sys.stderr)

    cdo = Cdo()

    # Apply mask and fldmean to each file individually before merging
    fldmean_files = []

    mask_file = "/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/mask_MEU-32.nc"

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

    variable_old = ""
    if variable == "tas":
        variable_old = "T_2M"
    elif variable == "tasmax":
        variable_old = "TMAX_2M"
    elif variable == "tasmin":
        variable_old = "TMIN_2M"
    elif variable == "pr":
        variable_old = "TOT_PREC"
    elif variable == "rsds":
        variable_old = "ASOD_S"
    elif variable == "sfcWind":
        variable_old = "SP_10M"
    else:
        raise ValueError("unknown variable")

    os.system(f"ncrename -v {variable_old},{variable} {dummy_data}")

    # At this point dummy_data is already masked+fldmean, so just do temporal averaging
    for resolution in temporal_resolution:
        if resolution == "yearly":
            cdo.yearmean(input=dummy_data, output=output_filename)
        elif resolution == "mon":
            cdo.monmean(input=dummy_data, output=output_filename)
        elif resolution == "day":
            cdo.daymean(input=dummy_data, output=output_filename)
        elif resolution == "1hr":
            os.system(f"cp {dummy_data} {output_filename}")
        else:
            raise ValueError("Unknown resolution or resolution not available.")


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
        for file in get_sorted_nc_files(folder):
            parts = os.path.basename(file).split("_")
            temp_resolution = os.path.dirname(file).split("/")[-1]

            # Navigate to the correct nested dict
            level = (
                info.setdefault(temp_resolution, {})
                .setdefault(parts[2], {})
                .setdefault(parts[3], {})
            )

            # Ensure the final key stores a list
            level.setdefault(parts[0], []).append(file)

    write_json_file(
        f"/work/bb1364/g260190_heinrich/UDAG/Data/json_files/{project}_{country}_{variable}_info.json",
        sort_dict_recursively(info),
    )


def precompute_masks(country):
    resolutions = {
        "MEU-3": "/work/bb1364/production_runs/work/IAEVALH01/post/yearly/T_2M/T_2M_1951010100-1951123123.ncz",
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
        mask_file = f"/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/mask_{res}2.nc"
        mask2d.to_netcdf(mask_file)

        saved_masks[res] = mask_file
        print(f"Saved mask for {res} to {mask_file}")

    return saved_masks

def sorted_resolution(list_of_wanted_resolutions) -> list[str]:
    temporary_resolutions = {
        "yearly": 1,
        "mon": 2,
        "day": 3,
        "1hr": 4,
    }
    
    return sorted(list_of_wanted_resolutions, key=lambda x: temporary_resolutions[x], reverse=True)

def main():
    variables = ["pr", "rsds", "sfcWind"]
    country = "Germany"
    project = "UDAG"

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

        for spatial_resolution in ["EUR-12"]:
            if (spatial_resolution == "CEU-3" and project != "NUKLEUS") or (
                spatial_resolution != "CEU-3" and project == "NUKLEUS"
            ):
                continue


            data_folders = find_folders(variable)
            for input_folder in data_folders:
                create_yearly_data(
                    input_folder,
                    output_folder,
                    overwrite,
                    list_of_wanted_resolutions,
                    variable,
                )

        create_info_json(output_folder)


if __name__ == "__main__":
    main()

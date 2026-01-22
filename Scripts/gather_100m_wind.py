import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import xarray as xr
from cdo import Cdo


def run_shell_command(command: str, time_minutes: int) -> None:
    try:
        subprocess.run(command, timeout=60 * time_minutes, shell=True)
    except subprocess.TimeoutExpired:
        print(f"The following command timed out: {command}", file=sys.stderr)


def calc_wind(mask_path, files)->list[str]:
    fld_means = []
    u_files = files

    v_files = [x.replace("ua100m", "va100m") for x in u_files]
    for u_wind, v_wind in zip(u_files, v_files):
        outfile = f"/scratch/g/g260190/sfcWind_100m_{hashlib.md5(u_wind.encode()).hexdigest()}"
        if os.path.exists(outfile):
            continue
        u_dummy = f"/scratch/g/g260190/u_{hashlib.md5(outfile.encode()).hexdigest()}.nc"
        v_dummy = f"/scratch/g/g260190/v_{hashlib.md5(outfile.encode()).hexdigest()}.nc"
        run_shell_command(f"touch {outfile}",1)
        run_shell_command(
            f"cdo -s -w -ifthen {mask_path} -sellonlatbox,0,8,54,57 {u_wind} {u_dummy}", 5
        )
        run_shell_command(
            f"cdo -s -w -ifthen {mask_path} -sellonlatbox,0,8,54,57 {v_wind} {v_dummy}", 5
        )
        run_shell_command(
            f"cdo -s -w -fldmean -expr,'sfcWind_100m=hypot(ua100m,va100m)' "
            f"-merge {u_dummy} {v_dummy} {outfile}",
            5,
        )
        fld_means.append(outfile)
    return fld_means


def write_json_file(filename: str, content: dict) -> None:
    with open(filename, "w") as file:
        json.dump(content, file, indent=4)


def get_highest_temporal_resolution(list_of_wanted_resolutions):
    temporary_resolutions = {
        "yearly": 1,
        "mon": 2,
        "day": 3,
        "1hr": 4,
    }

    return max(list_of_wanted_resolutions, key=lambda x: temporary_resolutions[x])


def find_folders(
    base, list_of_wanted_resolutions, variable, project, spatial_resolution
) -> list[str]:
    temporal_resolution = get_highest_temporal_resolution(list_of_wanted_resolutions)
    if temporal_resolution == "yearly":
        temporal_resolution = "mon"

    if project == "UDAG":
        base = Path(f"/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/{spatial_resolution}")
        pattern = f"CLMcom-*/*/*/*/ICON-CLM-202407-1-1/v1-r1/{temporal_resolution}/{variable}/v20240920"
    elif project == "NUKLEUS":
        base = Path("/work/bb1203/data_NUKLEUS_CMOR/CEU-3/")
        pattern = f"CLMcom-*/*/*/*/*/*/{temporal_resolution}/{variable}/*"
    else:
        raise ValueError(f"Unknown project: {project}")

    return sorted(str(m) for m in base.glob(pattern) if m.is_dir())


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
    if "ESGF_Buff" in folder:  # for UDAG
        return "_".join([parts[6], parts[7], parts[8], parts[9], variable]) + ".nc"
    elif "NUKLEUS" in folder:
        mapping = {
            "EC-Earth-Consortium-EC-Earth3-Veg": "EC-Earth3-Veg",
            "ECMWF-ERA5": "ERA5",
            "MIROC-MIROC6": "MIROC6",
            "MPI-M-MPI-ESM1-2-HR": "MPI-ESM1-2-HR",
        }
        if parts[6] not in mapping:
            raise ValueError("Cannot identify name of forcing")
        parts[6] = mapping[parts[6]]
        return "_".join([parts[4], parts[5], parts[6], parts[7], variable]) + ".nc"
    else:
        raise ValueError("Cannot identify project name")


def create_datasets(
    input_folder,
    output_folder_project,
    overwrite,
    list_of_wanted_resolutions,
    variable: str,
) -> None:
    highest_temporal_resolution = get_highest_temporal_resolution(
        list_of_wanted_resolutions
    )
    cdo = Cdo()
    dummy_data = ""
    for temporal_resolution in list_of_wanted_resolutions:
        output_folder = os.path.join(output_folder_project, temporal_resolution)
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        output_filename = os.path.join(
            output_folder, generate_filename(input_folder, variable)
        )
        if os.path.exists(output_filename) and not overwrite:
            continue

        if temporal_resolution == highest_temporal_resolution:
            files = get_sorted_nc_files(input_folder)
            dummy_data = "/scratch/g/g260190/dummy.nc"
            print(f"Working on {input_folder} {datetime.now()}", file=sys.stderr)

            # Apply mask and fldmean to each file individually before merging
            fldmean_files = []
            if "EUR-12" in input_folder:
                mask_file = "/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/North_Sea_EUR-12.nc"
            elif "MEU-3" in input_folder:
                mask_file = "/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/North_Sea_MEU-3.nc"
            else:
                raise ValueError(f"Unknown resolution in filename: {input_folder}")

            fldmean_files = calc_wind(mask_file, files)

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

        if not dummy_data:
            raise ValueError("input data is empty")
        # At this point dummy_data is already masked+fldmean, so just do temporal averaging
        if temporal_resolution == "yearly":
            cdo.yearmonmean(input=dummy_data, output=output_filename)
        elif temporal_resolution == "mon":
            cdo.monmean(input=dummy_data, output=output_filename)
        elif temporal_resolution == "day":
            cdo.daymean(input=dummy_data, output=output_filename)
        elif temporal_resolution == "1hr":
            run_shell_command(f"cp {dummy_data} {output_filename}",5)
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


def sorted_resolution(list_of_wanted_resolutions) -> list[str]:
    temporary_resolutions = {
        "yearly": 1,
        "mon": 2,
        "day": 3,
        "1hr": 4,
    }

    return sorted(
        list_of_wanted_resolutions, key=lambda x: temporary_resolutions[x], reverse=True
    )


def create_info_json(output_folder):
    print("Creating json files")
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


def precompute_masks():
    resolutions = {
        "EUR-12": "/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/EUR-12/CLMcom-Hereon/EC-Earth3-Veg/historical/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/fx/sftlf/v20240920/sftlf_EUR-12_EC-Earth3-Veg_historical_r1i1p1f1_CLMcom-Hereon_ICON-CLM-202407-1-1_v1-r1_fx.nc",
        "MEU-3": "/work/bb1149/ESGF_Buff/CORDEX-CMIP6/DD/MEU-3/CLMcom-Hereon/EC-Earth3-Veg/historical/r1i1p1f1/ICON-CLM-202407-1-1/v1-r1/fx/sftlf/v20240920/sftlf_MEU-3_EC-Earth3-Veg_historical_r1i1p1f1_CLMcom-Hereon_ICON-CLM-202407-1-1_v1-r1_fx.nc",
    }

    saved_masks = {}
    cdo = Cdo()
    for res, nc_file in resolutions.items():
        mask_file = (
            f"/work/bb1364/g260190_heinrich/UDAG/Scripts/masks/North_Sea_{res}.nc"
        )
        # os.system(f"cdo expr,'sftlf=(sftlf<50)?0:100' -sellonlatbox,0,8,54,57 {nc_file} {mask_file}")
        cdo.expr(
            "'sftlf=(sftlf<50)?0:100'",
            input=f"-sellonlatbox,0,8,54,57 {nc_file}",
            output=mask_file,
        )
        saved_masks[res] = mask_file
        print(f"Saved mask for {res} to {mask_file}")

    return saved_masks


def main():
    country = "North_Sea"
    project = "UDAG"

    list_of_wanted_resolutions = [
        "yearly",
        "mon",
        "day",
        "1hr",
    ]  # ["yearly", "mon", "day", "1hr"]
    list_of_wanted_resolutions = sorted_resolution(list_of_wanted_resolutions)

    overwrite = False
    precompute_masks()
    variable = "sfcWind_100m"
    output_folder = (
        f"/work/bb1364/g260190_heinrich/UDAG/Data/{project}/{country}/{variable}"
    )
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for spatial_resolution in ["EUR-12", "MEU-3", "CEU-3"]:
        if (spatial_resolution == "CEU-3" and project != "NUKLEUS") or (
            spatial_resolution != "CEU-3" and project == "NUKLEUS"
        ):
            continue
        data_folders = find_folders(
            project,
            list_of_wanted_resolutions,
            "ua100m",
            project,
            spatial_resolution,
        )
        for input_folder in data_folders:
            create_datasets(
                input_folder,
                output_folder,
                overwrite,
                list_of_wanted_resolutions,
                variable,
            )

    create_info_json(output_folder)


if __name__ == "__main__":
    main()

#%%
"""
Minimal BFP group-BrainSync atlas creation script.

This script intentionally stops after:
1. reading the config
2. selecting atlas subjects from the demographics CSV
3. loading/normalizing their BFP data
4. saving the atlas subject CSV files
5. running groupBrainSync
6. saving atlas.mat

Statistics, regression, hypothesis-testing, and visualization code are omitted.
"""

import os
import sys
import csv
import configparser

import scipy.io as spio


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

config_file = os.environ["CONFIG_FILE"]

config = configparser.ConfigParser()
config.read(config_file)

if not config.sections():
    raise RuntimeError(f"Could not read config file: {config_file}")


def get_config_value(name, default=None, required=True):
    """
    Find a config option by name regardless of which INI section contains it.
    This avoids loading the full statistics configuration just to run GroupSync.
    """
    for section in config.sections():
        if config.has_option(section, name):
            return config.get(section, name)

    if required:
        raise KeyError(
            f'Required config option "{name}" was not found in {config_file}'
        )

    return default


def as_bool(value):
    """Convert common config/CSV boolean representations to bool."""
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


# Only config fields required for atlas GroupSync
bfp_path = get_config_value("bfp_path")
out_dir = get_config_value("out_dir")
csv_fname = get_config_value("csv_fname")
data_dir = get_config_value("data_dir")
file_ext = get_config_value("file_ext")
colsubj = get_config_value("colsubj")
colvar_exclude = get_config_value("colvar_exclude")
colvar_atlas = get_config_value("colvar_atlas")
lentime = int(get_config_value("lentime"))
matcht = get_config_value("matcht")


# -------------------------------------------------------------------------
# Import only the BFP functions needed for atlas creation
# -------------------------------------------------------------------------

sys.path.append(os.path.join(bfp_path, "src", "stats"))
sys.path.append(os.path.join(bfp_path, "src", "BrainSync"))

from read_data_utils import load_bfp_dataT, write_text_timestamp
from brainsync import groupBrainSync


# -------------------------------------------------------------------------
# Output directory / log
# -------------------------------------------------------------------------

os.makedirs(out_dir, exist_ok=True)

log_fname = os.path.join(out_dir, "bfp_groupsync.log")

write_text_timestamp(
    log_fname,
    "Config file used: "
    + config_file
    + "\nAll outputs will be written in: "
    + out_dir,
)


# -------------------------------------------------------------------------
# Read demographics CSV and identify atlas subjects only
# -------------------------------------------------------------------------

subAtlas_IDs = []
subAtlas_fname = []
missing_files = []

with open(csv_fname, newline="") as csvfile:
    reader = csv.DictReader(csvfile)

    required_columns = {colsubj, colvar_exclude, colvar_atlas}
    missing_columns = required_columns.difference(reader.fieldnames or [])

    if missing_columns:
        raise KeyError(
            "Missing required demographics column(s): "
            + ", ".join(sorted(missing_columns))
        )

    for row in reader:
        subj = str(row[colsubj]).strip()

        if not subj:
            continue

        exclude = as_bool(row[colvar_exclude])
        use_for_atlas = as_bool(row[colvar_atlas])

        # GroupSync only needs non-excluded atlas subjects
        if exclude or not use_for_atlas:
            continue

        fname = os.path.join(
            data_dir,
            subj,
            "func",
            subj + file_ext,
        )

        if not os.path.isfile(fname):
            missing_files.append(fname)
            continue

        subAtlas_IDs.append(subj)
        subAtlas_fname.append(fname)


print(f"Atlas subjects selected: {len(subAtlas_IDs)}")

write_text_timestamp(
    log_fname,
    f"{len(subAtlas_IDs)} subjects will be used for Group BrainSync atlas creation.",
)

if missing_files:
    print(f"WARNING: {len(missing_files)} selected atlas files were not found.")
    write_text_timestamp(
        log_fname,
        f"WARNING: {len(missing_files)} selected atlas files were not found.",
    )

if len(subAtlas_IDs) < 2:
    raise RuntimeError(
        "Group BrainSync requires at least two usable atlas subjects."
    )


# -------------------------------------------------------------------------
# Load atlas-subject data
# -------------------------------------------------------------------------

print("Loading atlas subject data...")

subAtlas_data, subAtlas_numT = load_bfp_dataT(
    subAtlas_fname,
    lentime,
    matcht,
)


# -------------------------------------------------------------------------
# Save subjects used for atlas creation
# -------------------------------------------------------------------------

# Detailed table: participant ID + original number of timepoints
subjects_atlas_fname = os.path.join(out_dir, "subjects_atlas.csv")

with open(subjects_atlas_fname, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["participant_id", "numT"])
    writer.writerows(zip(subAtlas_IDs, subAtlas_numT))

# Simple ID-only list
subjects_atlas_ids_fname = os.path.join(out_dir, "subjects_atlas_IDs.csv")

with open(subjects_atlas_ids_fname, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["participant_id"])
    writer.writerows([[subj] for subj in subAtlas_IDs])

write_text_timestamp(
    log_fname,
    "Atlas subject list saved as "
    + subjects_atlas_fname
    + "\nAtlas subject ID list saved as "
    + subjects_atlas_ids_fname,
)


# -------------------------------------------------------------------------
# Run Group BrainSync
# -------------------------------------------------------------------------

write_text_timestamp(
    log_fname,
    "Running Group BrainSync algorithm for atlas creation.",
)

print("Running Group BrainSync...")

atlas_data, _, costdif, total_error = groupBrainSync(subAtlas_data)

print("Group BrainSync complete.")


# -------------------------------------------------------------------------
# Save atlas
# -------------------------------------------------------------------------

atlas_fname = os.path.join(out_dir, "atlas.mat")

spio.savemat(
    atlas_fname,
    {
        "atlas_data": atlas_data,
        "Costdif": costdif,
        "TotalError": total_error,
    },
)

write_text_timestamp(
    log_fname,
    "Atlas saved out as " + atlas_fname,
)

print(f"Atlas saved: {atlas_fname}")
print(f"Atlas subject table saved: {subjects_atlas_fname}")
print(f"Atlas subject ID list saved: {subjects_atlas_ids_fname}")

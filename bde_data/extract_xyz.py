"""Functions and a CLI to extract .xyz files from DES-typed csv data."""

# pylint: disable=too-many-arguments

import argparse
import os
from itertools import product
import pandas as pd
import numpy as np


def slice_des_data(data,
                   smiles0=None, smiles1=None, system_ids=None,
                   group_origs=None, group_ids=None, k_index=None,
                   geom_ids=None):
    """Slice a DES-type dataframe down according to various groupings

    Parameters
    ----------
    data : pandas DataFrame
        Data to be ingested as a pandas DataFrame or object with compatible
        methods. MUST contain the columns `elements` and `xyz`, as well as any
        columns utilized according to the optional parameters given (e.g.,
        `smiles0`, `group_id`, etc).
    smiles0, smiles1 : None or list of str, optional
        Lists of SMILES strings to be filtered for. If `smiles0` is not None
        but `smiles1` is None, then any dimers containing *a* monomer from
        `smiles0` will be retained. If `smiles0` is not None and `smiles1` is
        not None, then all dimers formed as the outer product of `smiles0` and
        `smiles1` will be
        retained.
    system_ids, group_ids, k_index, geom_ids: None or list of int, optional
        If not None, each list will be used to filter down entries according to
        the relevant column: `system_id`, `group_id`, and `geom_id`,
        respectively.  Only entries with a column matching an entry an entry in
        the relevant list will be retained.
    group_origs: None or list of str, optional
        If not None, then only entries with a `group_orig` matching an entry in
        `group_origs` will be retained.


    Returns
    -------
    data : pandas DataFrame
        A view of `data` that only contains entries matching one of the
        filters given in the args. This is view NOT a copy.
    """

    # Deal with easy selectors with sentinel values
    if system_ids is not None:
        data = data.loc[data.system_ids.isin(system_ids)]
    if group_origs is not None:
        data = data.loc[data.group_origs.isin(group_origs)]
    if group_ids is not None:
        data = data.loc[data.group_ids.isin(group_ids)]
    if k_index is not None:
        data = data.loc[data.k_index.isin(k_index)]
    if geom_ids is not None:
        data = data.loc[data.geom_ids.isin(geom_ids)]
    # SMILES inputs are special
    if smiles0 is None and smiles1 is None:  # Both have sentinel values
        pass
    elif smiles1 is None:  # Only 1 smiles set provided
        data = data.loc[data[['smiles0', 'smiles1']].isin(smiles0).any(axis=1)]
    else:  # Two smiles sets provided, time to get fancy
        # Reduce to rows-to-be-processed by doing the naive slice
        data = data.loc[data[['smiles0', 'smiles1']].isin(smiles0).any(axis=1)]
        data = data.loc[data[['smiles0', 'smiles1']].isin(smiles1).any(axis=1)]
        # We only want dimers that have one smiles from smiles0, and one smiles
        # from smiles1, but we may not know the order of them from the naive
        # product. Sorting gives us stable str comparisons.
        smiles = [tuple(sorted(x)) for x in product(smiles0, smiles1)]
        data = data.loc[data.apply(
            lambda x: tuple(sorted([x['smiles0'], x['smiles1']]))
            ).isin(smiles)]

    return data


def row_to_xyz(row):
    """Write a .xyz file based on one row of data, creating dirs if necessary

    The naming of an individual record will look like:
        smiles0_smiles1/group_orig/group_id/k_index.xyz

    The dirtree above the .xyz file will be created if required, using
    os.path.makedirs with exist_ok=True.

    The file will obey the standard .xyz format, with the comment line
    formatted as:
        smiles0, smiles1, group_orig, k_index, system_id, group_id, geom_id

    Parameters
    ----------
    row : pandas Series
        One row from a DES-type DataFrame. MUST contain the columns `elements`
        and `xyz`, as well as `smiles0`, `smiles1`, `system_id`, `group_orig`,
        `group_id`, and `k_index`.

    Returns
    -------
    None
    """
    xyz = np.fromstring(row.xyz, sep=' ', ).reshape(-1, 3)
    elements = row.elements.split()

    assert len(xyz) == len(elements), "Malformatted row! %s" % row

    lines = []
    lines.append(str(len(xyz)))  # n atoms
    lines.append(", ".join(  # comment line
        row[['smiles0', 'smiles1', 'group_orig', 'k_index', 'system_id',
             'group_id', 'geom_id']].astype(str)))
    for ixyz, ele in zip(xyz, elements):  # coordinates and atoms
        lines.append("%s %s %s %s" % (ele, *ixyz))

    # Write to-disk
    handle = os.path.join("%s_%s" % (row.smiles0, row.smiles1), row.group_orig,
                          str(row.group_id), "%i.xyz" % row.k_index)
    os.makedirs(os.path.dirname(handle), exist_ok=True)  # Create folder(s)
    with open(handle, "w") as stream:
        stream.write("\n".join(lines))


def extract_xyz(data, smiles0=None, smiles1=None, system_ids=None,
                group_origs=None, group_ids=None, k_index=None, geom_ids=None,
                verbose=False):
    """Create .xyz files from a dataframe with optional selecroot directory

    Given a DES-type dataframe, extract .xyz files for all entries matching
    various filtering criteria. The resulting files will be emplaced in a
    directory heirarchy that reflects the organization used in the DES data
    sets. The naming of an individual record will look like:
        smiles0_smiles1/group_orig/group_id/k_index.xyz

    The file will obey the standard .xyz format, with the comment line
    formatted as:
        smiles0, smiles1, group_orig, k_index, system_id, group_id, geom_id

    Parameters
    ----------
    data : pandas DataFrame
        Data to be ingested as a pandas DataFrame or object with compatible
        methods. MUST contain the columns `elements` and `xyz`, as well as any
        columns utilized according to the optional parameters given (e.g.,
        `smiles0`, `group_id`, etc).
    smiles0, smiles1 : None or list of str, optional
        Lists of SMILES strings to be filtered for. If `smiles0` is not None
        but `smiles1` is None, then any dimers containing *a* monomer from
        `smiles0` will be retained. If `smiles0` is not None and `smiles1` is
        not None, then all dimers formed as the outer product of `smiles0` and
        `smiles1` will be
        retained.
    system_ids, group_ids, k_index, geom_ids: None or list of int, optional
        If not None, each list will be used to filter down entries according to
        the relevant column: `system_id`, `group_id`, and `geom_id`,
        respectively.  Only entries with a column matching an entry an entry in
        the relevant list will be retained.
    group_origs: None or list of str, optional
        If not None, then only entries with a `group_orig` matching an entry in
        `group_origs` will be retained.
    verbose: bool, optional
        Toggle to print a running count of records found and processed


    Returns
    -------
    None

    """
    data_view = slice_des_data(data, smiles0, smiles1, system_ids, group_origs,
                               group_ids, k_index, geom_ids)
    for i, row in data_view.iterrows():
        if verbose:
            print("Processing record %i/%i" % (i, len(data_view)), end='\r')
        row_to_xyz(row)


def main():
    """Provide a CLI interface"""
    parser = argparse.ArgumentParser(
        description="""Produce .xyz files from a DES-typed csv.

Files will be located in a dirtree that looks like:
    args.output/smiles0_smiles1/group_orig/group_id/k_index.xyz
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("csv", type=os.path.abspath,
                        help="DES-typed CSV, e.g., DES370K.csv")
    parser.add_argument("output", type=os.path.abspath,
                        help="Output directory for created dirtree. "
                             "Will be created if it doesn't exist")
    parser.add_argument("--smiles0", type=str, nargs='+', default=None,
                        help="SMILES strings to filter on. If used, only "
                             "dimers containing *a* monomer in `smiles0` will "
                             "be kept.")
    parser.add_argument("--smiles1", type=str, nargs='+', default=None,
                        help="SMILES strings to filter on. If used, only "
                             "dimers matching the product of `smiles0` and "
                             "`smiles1` will be kept. Ignored if `smiles0` is "
                             "not set.")
    parser.add_argument("--system_ids", type=int, nargs='+', default=None,
                        help="system_id(s) to filter on.")
    parser.add_argument("--group_origs", type=str, nargs='+', default=None,
                        choices=['qm_opt_dimer', 'md_solvation', 'md_dimer',
                                 'md_nmer'],
                        help="group_orig(s) to filter on.")
    parser.add_argument("--group_ids", type=int, nargs='+', default=None,
                        help="group_id(s) to filter on.")
    parser.add_argument("--k_index", type=int, nargs='+', default=None,
                        help="group_index(s) to filter on.")
    parser.add_argument("--geom_ids", type=int, nargs='+', default=None,
                        help="geom_id(s) to filter on.")
    parser.add_argument("--verbose", action='store_true',
                        help="Show a running counter of rows processed")

    args = parser.parse_args()

    # Set up our root directory
    os.makedirs(args.output, exist_ok=True)
    cwd = os.path.abspath(os.curdir)
    os.chdir(args.output)

    data = pd.read_csv(args.csv)

    extract_xyz(data, args.smiles0, args.smiles1, args.system_ids,
                args.group_origs, args.group_ids, args.k_index, args.geom_ids,
                args.verbose)

    os.chdir(cwd)


if __name__ == "__main__":
    main()

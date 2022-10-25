#!/opt/conda/envs/sc-dandelion-container/bin/python
import argparse
import dandelion as ddl
import os


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5ddl",
        required=True,
        help=("Dandelion object to call changeo clonotypes on."),
    )
    parser.add_argument(
        "--manual_threshold",
        help=("Optional manual override of SHazaM threshold."),
    )
    parser.add_argument(
        "--plot_file",
        help=(
            "File name to save PDF of SHazaM plot to. Defaults to the input object name "
            + "with _shazam.pdf appended."
        ),
    )
    parser.add_argument(
        "--key_added",
        default="changeo_clone_id",
        help=(
            "The column name to store the identified clone ID under in the object. "
            + "Defaults to changeo_clone_id."
        ),
    )
    parser.add_argument(
        "--h5ddl_out",
        help=(
            "Path to save Dandelion object with changeo clonotypes to. Defaults to the "
            + "input object name with _changeo.h5ddl appended."
        ),
    )
    args = parser.parse_args()
    # set up the default names of files if need be. needs the base name of the file
    basename = os.path.splitext(args.h5ddl)[0]
    if args.plot_file is None:
        args.plot_file = basename + "_shazam.pdf"
    if args.h5ddl_out is None:
        args.h5ddl_out = basename + "_changeo.h5ddl"
    return args


def main():
    # parse arguments
    args = parse_args()
    # the actual process is easy. the dependencies quite a bit less so
    vdj = ddl.read_h5ddl(args.h5ddl)
    ddl.pp.calculate_threshold(
        vdj,
        manual_threshold=float(args.manual_threshold),
        save_plot=args.plot_file,
    )
    ddl.tl.define_clones(vdj, key_added=args.key_added)
    vdj.write(args.h5ddl_out)


if __name__ == "__main__":
    main()

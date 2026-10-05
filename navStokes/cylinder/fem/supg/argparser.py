import argparse 

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-y",
        "--yaml",
        type=str,
        default="config.yaml",
        help="Name of configuration file, in .yaml format."
    )
    parser.add_argument(
        "-m",
        "--mesh",
        type=str,
        default="quad2D_2504cells.msh",
        help="Name of meshfile to be read, resolved against cylinder/mesh/ -- "
        "see cylinder/mesh/meshes_summary.csv for what's available."
    )
    parser.add_argument(
        "-c",
        "--case",
        type=int,
        default=3,
        choices=[1, 2, 3],
        help="DFG case number -- see the 'cases' section of config.yaml "
        "for each one's name/parameters.",
    )
    return parser.parse_args()


#########################################################################################################################
# usage: generate_trajectories.py [-h] -s SNAPSHOTS -o OUTPUT [-e FILE_EXT] [--radial [RADIAL]] [--kappa [KAPPA]]       #
#                                                                                                                       #
# Generate Node trajectories from snapshots' coordinates files                                                          #
#                                                                                                                       #
# optional arguments:                                                                                                   #
#   -h, --help            show this help message and exit                                                               #
#   -s SNAPSHOTS, --snapshots SNAPSHOTS                                                                                 #
#                         The folder with the network snapshot coordinates in the form <string_snapID.inf_coord>        #
#   -o OUTPUT, --output OUTPUT                                                                                          #
#                         The output folder to store the node trajectories. Output files names are the node Ids         #
#                         containing space separated lines with the following fields <snapshot_id, value>               #
#   -e FILE_EXT, --file_ext FILE_EXT                                                                                    #
#                         The file extension of the network snapshots input files                                       #
#   --radial [RADIAL]     To extract radial instead of angular coordinates. Remove this option to extract the angular   #
#                         coordinates of the nodes.                                                                     #
#   --kappa [KAPPA]       To extract kappa instead of angular coordinates. Remove this option to extract the angular or #
#                         radial coordinates of the nodes.                                                              #
#                                                                                                                       #
# Authors: Costas Iordanou                                                                                              #
#########################################################################################################################

import argparse
import os
import re


# To print custom exceptions
class tool_exception(Exception):
    pass


# Extract integer number from string
def extract_numbers_from_string(str):
    numList = re.findall(r'\d+', str)
    # We do not accept file names with more than one integer number in 
    # their filename since we are unable to identify the snapshot number
    if len(numList) > 1:
        raise tool_exception("{} file name has more than one number".format(str))
    return numList[0]


# Get the total number of snapshots from file names for a given folder
def get_number_of_snapshots(path, f_ext=".inf_coord"):
    snaps = []
    for f in os.listdir(path):
        if f.endswith(f_ext):
            snap_id = int(extract_numbers_from_string(f))
            # Here we record the snapshot id and the filename in 
            # order to be able to detect gaps between snapshots
            # i.e. snapshots that do not increase by 1
            snaps.append({"id": snap_id, "f_name": f})
    return snaps


# To handle boolean arguments
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Read Mercator rotated output files
def load_raw_mercator(f_name, radial, kappa):
    nodes = {}
    with open(f_name, "r") as fp:
        for l in fp.readlines():
            # Clean up the line from tabs and newline
            line = re.sub(' +', ' ', l.strip())
            # Ignore comment lines
            if not line.startswith("#"):
                # Split line to the different fields:
                #        Vertex = id       Inf.Kappa       Inf.Theta = theta_rad    Inf.Hyp.Rad.
                tl = line.split(" ")
                if len(tl) != 4:
                    raise tool_exception("Input file fields mismatch. Expected 4 found {}. Make sure you are using the correct network type.".format(len(tl)))
                # Grab the node id and the associated coordinate
                if tl[0] not in nodes:
                    if radial:
                        nodes[tl[0]] = float(tl[3])
                    elif kappa:
                        nodes[tl[0]] = float(tl[1])
                    else:
                        nodes[tl[0]] = float(tl[2])
                else: # Sanity check just to make sure that all went as expected
                    if radial:
                        print("Node {} has already radial: {} - New: {}".format(tl[0], nodes[tl[0]], tl[3]))
                    elif kappa:
                        print("Node {} has already kappa: {} - New: {}".format(tl[0], nodes[tl[0]], tl[1]))
                    else:
                        print("Node {} has already theta: {} - New: {}".format(tl[0], nodes[tl[0]], tl[2]))
                    raise tool_exception("Same node on the same snapshot. Something went wrong. Exiting...")

    return nodes


# Extract node trajectories from snapshots
def extract_nodes_trajectories(snaps):
    # Get the ordered (snapshot) node coordinates 
    node_trajectories = {}
    for s in snaps:
        nodes = s['nodes']
        for n in nodes:
            val = float(nodes[n])
            if n not in node_trajectories:
                node_trajectories[n] = []
            node_trajectories[n].append({s['id']: val})
    return node_trajectories


# Write final data to files
def write_trajectories_files(out_folder, nodes_trajectories):
    
    # Check if output folder exists
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    f_count = 0
    total_files = len(nodes_trajectories)
    for n in nodes_trajectories:
        lines = []
        lines.append("#snapId value")
        # Sort node snapshots
        ord_snaps = sorted(nodes_trajectories[n], key=lambda d: list(d.keys()))
        for e in ord_snaps:
            snap = list(e.keys())[0]
            val = list(e.values())[0]
            lines.append("{} {}".format(snap, val))
        f_name = os.path.join(out_folder, "Node_{}.txt".format(n))
        try:
            with open(f_name, "w", newline="") as fp:
                fp.writelines("\n".join(lines))
        except Exception as ex:
            print(f"Fail to save file: {f_name}\nException: {ex}")
        print("Saving node {} out of {}".format(f_count, total_files), end="\r", flush=True)
        f_count += 1
    print("Saving node {} out of {}".format(f_count, total_files))


if __name__ == "__main__":
    # Handle the script arguments
    parser = argparse.ArgumentParser(description="Generate Node trajectories from snapshots' coordinates files", epilog="Authors: Costas Iordanou")
    parser.add_argument('-s', '--snapshots', 
        help='The folder with the network snapshot coordinates in the form <string_snapID.inf_coord>', 
        required=True
    )
    parser.add_argument('-o', '--output', 
        help='The output folder to store the node trajectories. Output files names are the node Ids containing space separated lines with the following fields <snapshot_id, value>', 
        required=True
    )
    parser.add_argument('-e', '--file_ext', help="The file extension of the network snapshots input files", default='.inf_coord', required=False)
    parser.add_argument('--radial', type=str2bool, const=True, nargs='?', default=False, help="To extract radial instead of angular coordinates. \
        Remove this option to extract the angular coordinates of the nodes.")
    parser.add_argument('--kappa', type=str2bool, const=True, nargs='?', default=False, help="To extract kappa instead of angular coordinates. \
        Remove this option to extract the angular or radial coordinates of the nodes.")

    args =parser.parse_args()
    input_folder = args.snapshots
    output_folder = args.output
    f_ext = args.file_ext
    radial = args.radial
    kappa = args.kappa
    
    if radial and kappa:
        raise tool_exception("You need to select either Radial or kappa coordinates but not both. Exiting...")

    if input_folder == output_folder:
        raise tool_exception("Input folder is equal to the output folder. Please select a different folder.\nExiting...")

    ### Read network snapshots ###
    snaps = get_number_of_snapshots(input_folder, f_ext=f_ext)
    total_snaps = len(snaps)

    # Check if we have valid snapshots
    if total_snaps == 0:
        raise tool_exception("No snapshots found. Make sure you point to a folder with network snapshot files (.inf_coord).\nExiting...")

    ### Read snapshots input files ###
    print("Reading snapshots input files...")
    f_count = 0
    for s in snaps:
        print("Reading file: {} out of {}".format(f_count, total_snaps), end="\r", flush=True)
        f_name = os.path.join(input_folder, s['f_name'])
        s['nodes'] = load_raw_mercator(f_name=f_name, radial=radial, kappa=kappa)
        f_count += 1
    print("Reading file: {} out of {}".format(f_count, total_snaps))

    ### Extract node trajectories ###
    print("Extracting nodes trajectories...")
    nodes_trajectories = extract_nodes_trajectories(snaps)
    total_nodes = len(nodes_trajectories)
    print("Total nodes found: {}".format(total_nodes))

    ### Write results to files ###
    write_trajectories_files(out_folder=output_folder, nodes_trajectories=nodes_trajectories)

    ### Print final report ###
    print("\nNodes Trajectory Report\n{}\nInput folder: {}\nOutput folder: {}\nTotal processed snapshots: {}\nTotal Nodes: {}\n{}\n".format(
        "="*60,input_folder, output_folder, total_snaps, total_nodes, "="*60)
    )
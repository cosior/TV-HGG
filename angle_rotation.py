########################################################################################################################
# usage: angle_rotation.py [-h] -s SNAPSHOTS -o OUTPUT                                                                 #
#                                                                                                                      #
# Rotate Mercator angular coordinates                                                                                  #
#                                                                                                                      #
# optional arguments:                                                                                                  #
#  -h, --help            show this help message and exit                                                               #
#  -s SNAPSHOTS, --snapshots SNAPSHOTS                                                                                 #
#                        The folder with the Mercator snapshots in the form <string_snapID.inf_coord>                  #
#  -o OUTPUT, --output OUTPUT                                                                                          #
#                        The output folder to store the rotated angular coordinates. Output files are space separated  #
#                        files with the following fields <id Kappa Theta Hyp.Rad> similar to the Mercator output files #
#                                                                                                                      #
# Authors: Costas Iordanou                                                                                             #
########################################################################################################################

import argparse
import textwrap
import os
import math
import re

PI = math.pi


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


# Read Mercator output files
def load_raw_mercator(f_name):
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
                # Grab the node id and the associated coordinates
                node = tl[0]
                try:
                    node = int(tl[0])
                except:
                    pass
                if node not in nodes:
                    nodes[node] = (float(tl[1]), float(tl[2]), float(tl[3]))
                else: # Sanity check just to make sure that all went as expected
                    print("Same node {} on the same snapshot. Something went wrong. Exiting...".format(tl))
                    exit()
    return nodes

# The main rotation function
def rotate_snapshots(first_file, second_file, output_file):
    
    # Read the first file
    f_nodes = load_raw_mercator(first_file)
    x_1 = {}
    y_1 = {}
    k = {}
    for n in f_nodes:
        x_1[n] = math.cos(f_nodes[n][1])
        y_1[n] = math.sin(f_nodes[n][1])
        k[n] = f_nodes[n][0]

    # Read the second file
    s_nodes = load_raw_mercator(second_file)
    w_1 = {}
    z_1 = {}
    w_1_t = {}
    z_1_t = {}
    r = {}
    tildek = {}
    for n in s_nodes:
        r[n] = s_nodes[n][2]
        kagg = s_nodes[n][0]
        theta = s_nodes[n][1]
        w_1[n] = math.cos(theta)
        z_1[n] = math.sin(theta)
        tildek[n] = kagg
        # Check also with 2*p1 translation first.
        theta = 2*PI - s_nodes[n][1]
        w_1_t[n] = math.cos(theta)
        z_1_t[n] = math.sin(theta)

    ### Find angle of rotation ###
    # without translation
    sum_1 = 0
    sum_2 = 0
    for n in sorted(w_1.keys()):
        if n in x_1:
            sum_1 += w_1[n]*y_1[n] - z_1[n]*x_1[n]
            sum_2 += w_1[n]*x_1[n] + z_1[n]*y_1[n]
    r_theta = math.atan2(sum_1, sum_2)

    # With translation first
    sum_1 = 0
    sum_2 = 0
    for n in sorted(w_1_t.keys()):
        if n in x_1:
            sum_1 += w_1_t[n]*y_1[n] - z_1_t[n]*x_1[n]
            sum_2 += w_1_t[n]*x_1[n] + z_1_t[n]*y_1[n]
    r_theta_t = math.atan2(sum_1, sum_2)

    diff = 0
    mse = 0
    for n in sorted(w_1.keys()):
        if n in x_1:
            u_1 = w_1[n]*math.cos(r_theta) - z_1[n]*math.sin(r_theta)
            v_1 = w_1[n]*math.sin(r_theta) + z_1[n]*math.cos(r_theta)

            r_theta_1 = math.atan2(y_1[n], x_1[n])
            if r_theta_1 < 0:
                r_theta_1 = 2*PI + r_theta_1
            r_theta_2 = math.atan2(v_1, u_1)
            if r_theta_2 < 0:
                r_theta_2 = 2*PI + r_theta_2
            diff += abs(r_theta_1-r_theta_2)
            mse += (r_theta_1-r_theta_2)**2

    diff_t = 0
    mse_t = 0
    for n in sorted(w_1_t.keys()):
        if n in x_1:
            u_1 = w_1_t[n]*math.cos(r_theta_t) - z_1_t[n]*math.sin(r_theta_t)
            v_1 = w_1_t[n]*math.sin(r_theta_t) + z_1_t[n]*math.cos(r_theta_t)

            r_theta_1_t = math.atan2(y_1[n], x_1[n])
            if r_theta_1_t < 0:
                r_theta_1_t = 2*PI + r_theta_1_t
            r_theta_2_t = math.atan2(v_1, u_1)
            if r_theta_2_t < 0:
                r_theta_2_t = 2*PI + r_theta_2_t
            diff_t += abs(r_theta_1_t-r_theta_2_t)
            mse_t += (r_theta_1_t-r_theta_2_t)**2

    LINES = []
    if diff < diff_t:
        for n in sorted(w_1.keys()):
            u_1 = w_1[n]*math.cos(r_theta) - z_1[n]*math.sin(r_theta)
            v_1 = w_1[n]*math.sin(r_theta) + z_1[n]*math.cos(r_theta)
            
            r_theta_2 = math.atan2(v_1, u_1)
            if r_theta_2 < 0:
                r_theta_2 = 2*PI + r_theta_2
            LINES.append("{} {:5f} {:5f} {:5f}".format(n, tildek[n], r_theta_2, r[n]))
    else:
        for n in sorted(w_1_t.keys()):
            u_1 = w_1_t[n]*math.cos(r_theta_t) - z_1_t[n]*math.sin(r_theta_t)
            v_1 = w_1_t[n]*math.sin(r_theta_t) + z_1_t[n]*math.cos(r_theta_t)
            r_theta_2_t = math.atan2(v_1, u_1)
            if r_theta_2_t < 0:
                r_theta_2_t = 2*PI + r_theta_2_t
            LINES.append("{} {:5f} {:5f} {:5f}".format(n, tildek[n], r_theta_2_t, r[n]))

    with open(output_file, "w") as fp:
        fp.writelines("\n".join(LINES))


# The main snapshots loop function
def main_loop(input_folder, output_folder):
    
    # Get the file names of a given folder
    files = get_number_of_snapshots(input_folder)
    
    # Sort the files based on their snapshot number
    sorted_files = sorted(files, key=lambda d: d['id'])
    
    # Get the current working directory
    cwd = os.getcwd()
    
    # Rotate the first snapshot for format consistency
    # The coordinates remain the same but the format will be the same as the rest of the rotated files
    current_f = sorted_files[0] # The initial network snapshot
    next_f = sorted_files[0] # Again here we need the initial network snapshot
    ref_coords = os.path.join(cwd, input_folder, current_f["f_name"]) # Build the reference coordinates file path
    toRotate = os.path.join(cwd, input_folder, next_f['f_name']) # Build the file path for the snapshot to be rotated (initially is just the reference snapshot)
    out_f = os.path.join(cwd, output_folder, "r_{}.inf_coord".format(next_f["id"])) # Build the output file path
    
    # Inform the user and execute the rotation
    print("Rotating snapshot {} to {}".format(next_f['f_name'], out_f.split("\\")[-1]))
    rotate_snapshots(ref_coords, toRotate, out_f)

    # The output of this step is going to be the reference snapshot for the next one.
    next_ref = out_f
    
    # The main rotation loop
    for i in range(len(sorted_files)-1):
        current_f = sorted_files[i] 
        next_f = sorted_files[i+1]
        
        ref_coords = next_ref # The file path of the reference snapshot (already rotated)
        toRotate = os.path.join(cwd, input_folder, next_f['f_name']) # The file path of the next snapshot to be rotated 
        out_f = os.path.join(cwd, output_folder, "r_{}.inf_coord".format(next_f["id"])) # The file path of the output file

        # Inform the user and execute the rotation
        print("Rotating snapshot {} to {}".format(next_f['f_name'], out_f.split("\\")[-1]))
        rotate_snapshots(ref_coords, toRotate, out_f)

        # The output of this step is going to be the reference snapshot for the next one.
        next_ref = out_f


if __name__ == "__main__":
    # Handle the script arguments
    parser = argparse.ArgumentParser(description="Rotate Mercator angular coordinates",  formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''    
      Authors: Costas Iordanou
      '''))
    parser.add_argument('-s', '--snapshots', 
        help='The folder with the Mercator snapshots in the form <string_snapID.inf_coord>', 
        required=True
    )
    parser.add_argument('-o', '--output', 
        help='The output folder to store the rotated angular coordinates. Output files are space separated files with the following fields <id Kappa Theta Hyp.Rad> similar to the Mercator output files', 
        required=True
    )
    args =parser.parse_args()
    input_folder = args.snapshots
    output_folder = args.output

    # Just make sure we do not accidentally mix the rotated with the un-rotated coordinates
    if input_folder == output_folder:
        raise tool_exception("Input folder is equal to the output folder. This will overwrite original coordinate files.\nExiting...")

    # Check if the output folder exists and if not create it.
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # The main rotation function
    main_loop(input_folder=input_folder, output_folder=output_folder)

    # Check if all input files are actually available in rotated format
    in_files = len(get_number_of_snapshots(input_folder))
    out_files = len(get_number_of_snapshots(output_folder))
    if in_files != out_files:
        print("We have {} files that are not rotated.".format(in_files-out_files))

    # Show final report to the user
    print("\nRotation Report\n{}\nTotal rotated snapshot: {}\nOutput path: {}\n{}\n".format(
        "="*60, len(get_number_of_snapshots(output_folder)), os.path.join(os.getcwd(), output_folder), "="*60
    ))

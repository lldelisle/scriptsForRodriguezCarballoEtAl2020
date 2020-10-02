import cooler
import os
import pandas as pd
import argparse
import sys

def build_cool(input_bedgraph, frag_per_bin, pixels_file, output_file):
    # Build a data frame from the input_bedgraph
    input_df = pd.read_csv(input_bedgraph, sep='\t', header = None)
    # Using the input bedgraph and the number of fragment en each bead
    # We can deduce the coordinates of each bead
    data = {'chrom':input_df[0][0::frag_per_bin].tolist(),
            'start':input_df[1][0::frag_per_bin].tolist(),
            'end':input_df[1][frag_per_bin::frag_per_bin].tolist()}

    # We may need to adjust the final coordinate
    if len(data['end']) < len(data['start']):
        data['end'].append(input_df[2][len(input_df[2]) - 1])

    # We store the data in a panda dataframe
    new_df = pd.DataFrame(data, columns=['chrom', 'start', 'end'])

    # We store in panda dataframe
    pixels = pd.read_csv(pixels_file, header = None)
    pixels.columns = ['bin1_id', 'bin2_id', 'count']

    # We convert count to integer as this is the regular format of cool file
    pixels['count'] = pixels['count'].astype(int)

    # We only keep upper corner
    pixels = pixels[pixels.bin1_id <= pixels.bin2_id]

    # We need to add a pixel starting at 0 to be 'real' cool file:
    new_coos = pd.DataFrame({'chrom':[data['chrom'][0]],
                             'start':[0],
                             'end':[data['start'][0]]},
                             columns=['chrom', 'start', 'end']).append(new_df, ignore_index=True)

    # We need to update the bins_id in the pixels:
    pixels['bin1_id'] +=1
    pixels['bin2_id'] +=1

    # We create the cool file
    cooler.create_cooler(output_file, new_coos, pixels)


argp = argparse.ArgumentParser(
    description=("Convert virtual HiC file (pixel_file)"
                 " in cool file."))
argp.add_argument('--data_dir', default=None,
                  help="data_dir used in 4Cin where input bedgraph are put.")
argp.add_argument('--working_dir', default=None,
                  help="working_dir used in 4Cin where you have the result.")
argp.add_argument('--prefix', default="HoxD",
                  help="prefix used in 4Cin, if not provided, the name of the pixels_file should be vhic_'prefix'.txt.")
argp.add_argument('--input_bedgraph', default=None,
                  help="A bedgraph which has been used as input to 4Cin.")
argp.add_argument('--frag_per_bin', default=None, type=int,
                  help="The number of fragment per bin (this information is in the log.txt).")
argp.add_argument('--pixel_file', default=None,
                  help="The pixel file to transform to cool (usually vhic_'prefix'.txt).")
argp.add_argument('--output_file', default=None, required=True,
                  help="The output file cool.")

args = argp.parse_args()

if os.path.exists(args.output_file):
  raise Exception("The output file already exists. Please remove it or choose another output file.")

# We need a bedgraph for the coordinates
if args.input_bedgraph is None:
    # We need to guess it
    if args.data_dir is None:
        raise Exception("input_bedgraph is not given nor data_dir, cannot guess the coordinates.")
    else:
        # We get one input bedgraph to get the coordinates of the beads in 4Cin model
        # We assume the directory contains the input bedgraphs,
        # maybe the working_dir,
        # primers.txt and maybe a directory named forigv but nothing else
        if not os.path.exists(args.data_dir):
            raise Exception("data_dir provided does not exists")
        data_dir = args.data_dir
        if data_dir[-1] == "/":
            data_dir = data_dir[:-1]
        if args.working_dir is not None and args.working_dir.startswith(data_dir):
            relative_working_dir = args.working_dir.replace(data_dir + '/', '').split('/')[0]
        else:
            relative_working_dir = 'wd'
        input_bedgraph_relative = [f for f in os.listdir(data_dir) if f not in [relative_working_dir, 'primers.txt', 'forigv']][0]
        input_bedgraph = os.path.join(data_dir, input_bedgraph_relative)
else:
    input_bedgraph = args.input_bedgraph

# We need the number of fragments per bead:
if args.frag_per_bin is None:
    # We need to find it in the log file
    if args.working_dir is None:
        print("The working_dir is not specified, we assume it is part of data_dir")
        if args.data_dir is None:
            raise Exception("frag_per_bin is not given nor working_dir nor data_dir...")
        directory_to_use = args.data_dir
        log_file = None
    else:
        directory_to_use = args.working_dir
        if args.prefix is None:
            log_file = None
        else:
            log_file = os.path.join(args.working_dir, args.prefix, 'log.txt')
            if not os.path.exists(log_file):
                print("The log.txt was not found in " + log_file + " will try to find it in the working_dir")
                log_file = None
    if log_file is None:
        log_files = [os.path.join(path, 'log.txt') for (path, dirs, files) in os.walk(directory_to_use) if 'log.txt' in files]
        if len(log_files) == 0:
            raise Exception("frag_per_bin is not given and cannot find a log.txt in the working_dir.")
        if len(log_files) > 1:
            raise Exception("frag_per_bin is not given and found multiple log.txt in the working_dir.")
        log_file = log_files[0]

    frag_per_bin = None
    # In the log.txt file there are the number of fragments in each bead
    with open(log_file, 'r') as f:
        for line in f:
            if 'Fragments in each bead' in line:
                frag_per_bin = int(line.strip().split(':')[1])
    
    if frag_per_bin is None:
        raise Exception("frag_per_bin is not given and this information was not in the log file.")
else:
    frag_per_bin = args.frag_per_bin

# We need the pixel file:
if args.pixel_file is None:
    # We need to find it.
    if args.working_dir is None:
        print("The working_dir is not specified, we assume it is part of data_dir")
        if args.data_dir is None:
            raise Exception("pixel_file is not given nor working_dir nor data_dir...")
        directory_to_use = args.data_dir
        pixel_files = None
    else:
        directory_to_use = args.working_dir
        if args.prefix is None:
            pixel_files = None
        else:
            # The name of the final directory is not fixed but contains final
            relative_final_dirs = [ d for d in os.listdir(os.path.join(args.working_dir, args.prefix)) if 'final' in d ]
            potential_pixel_files = [os.path.join(args.working_dir, args.prefix, final_dir, 'vhic_' + args.prefix + '.txt') for final_dir in relative_final_dirs]
            pixel_files = [pixel_file for pixel_file in potential_pixel_files if os.path.exists(pixel_file)]
            if len(pixel_files) == 0:
                print("No pixel_file was not found in " + os.path.join(args.working_dir, args.prefix) + " will try to find it in the working_dir")
                pixel_files = None
    if pixel_files is None:
        pixel_files = [os.path.join(path, f) for (path, dirs, files) in os.walk(directory_to_use) for f in files if f.startswith('vhic_') and f.endswith('.txt')]
    if len(pixel_files) == 0:
        raise Exception("pixel_file is not given and cannot find a vhic_*.txt in the working_dir.")
    if len(pixel_files) > 1:
        raise Exception("pixel_file is not given and found multiple vhic_*.txt in the working_dir.")
    pixel_file = pixel_files[0]
else:
    pixel_file = args.pixel_file
    
build_cool(input_bedgraph, frag_per_bin, pixel_file, args.output_file)




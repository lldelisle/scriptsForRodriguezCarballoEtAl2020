import sys
import os

import cooler
import scipy.stats
import numpy as np
import pyBigWig
import gzip

path_GitHub_table = sys.argv[1]
path_4Cin = sys.argv[2]
path_4C_geo = sys.argv[3]
bed_file_to_quantify = sys.argv[4]
table_to_quantify = sys.argv[5]
output_file_stats_4C = sys.argv[6]
output_file_quantif_4C = sys.argv[7]

plotted_region = (74401941, 75800320)

subTADs_file = os.path.join(path_GitHub_table, "subTADs.bed")

# First I load the subTADs coordinates from bed
subTADs={}
subTADs_names = []

with open(subTADs_file, 'r') as f:
    for line in f:
        if line[:5] == 'track':
            continue
        ls = line.strip().split("\t")
        subTADs[ls[3]] = f"{ls[0]}:{ls[1]}-{ls[2]}"
        subTADs_names.append(ls[3])

# Then I write the comparisons to evaluate in the different figures
pairs = {'fig2': ('E9_FLB_delCS3840', 'E9_FLB_wt_like_delCS3840'),
         'figS4': ('E12_PFL_delCTCFs', 'E12_PFL_wt'),
         'fig3': ('E9_FLB_invCS3840', 'E9_FLB_wt_like_invCS3840')}

# For each figure
for fig in pairs.keys():
    # The correct cool matrices are loaded:
    mutant_cool = path_4Cin + "/average__" + pairs[fig][0] + "__all.cool"
    wt_cool = path_4Cin + "/average__" + pairs[fig][1] + "__all.cool"

    both_cool_files = [wt_cool, mutant_cool]
    both_cool = [cooler.Cooler(f) for f in both_cool_files]

    # For each comparison: subTAD1 vs subTAD1, subTAD1 vs subTAD2, subTAD2 vs subTAD2
    for i in range(len(subTADs_names)):
        for j in range(i, len(subTADs_names)):
            print("\n" + fig)
            print(f"{subTADs_names[i]} vs {subTADs_names[j]}")
            # Data is loaded from cool file
            data = [c.matrix(balance=False).fetch(subTADs[subTADs_names[i]], subTADs[subTADs_names[j]]) for c in both_cool]
            # When i and j are identical only the upper matrix is kept
            if i == j:
                data = [d[np.triu_indices(d.shape[0])] for d in data]
            else:
                data = [d.flatten() for d in data]
            # Mean distance is plotted
            print([np.mean(d) for d in data])
            # The fold-change in mean
            print("Fold-change of mean: {}".format(np.mean(data[1])/np.mean(data[0])))
            # U test is performed
            print(scipy.stats.mannwhitneyu(*data))

# Stats on 4C: nb of fragments covered, signal in trans etc.
with open(output_file_stats_4C, 'w') as fo:
    fo.write("file_name\tNb_frags\tTrans\tSignalInPlottedRegion\tFragmentsInPlottedRegion\n")
    # The segToFrag contains the scores stored as bw file normalized to million mapped reads
    for segToFrag in [f for f in os.listdir(path_4C_geo) if f.startswith("segToFrag") and f.endswith(".bw")]:
        # Open the bigwig
        bw = pyBigWig.open(os.path.join(path_4C_geo, segToFrag))
        all_count = 0
        all_signal = 0
        signal_chr2 = 0
        # For each chromosome we sum the signal and count the number of fragments
        for chrom in bw.chroms():
            intervals = bw.intervals(chrom, 0, bw.chroms(chrom))
            all_count += len(intervals)
            signal = sum([v[2] for v in intervals])
            all_signal += signal
            if chrom == "chr2":
                signal_chr2 = signal
                signal_plotted_region = sum([v[2] for v in bw.intervals(chrom, *plotted_region)])
                frags_plotted_region = len(bw.intervals(chrom, *plotted_region))
        fo.write(f"{segToFrag}\t{all_count}\t{100 * (1 - signal_chr2 / all_signal)}\t{100 * signal_plotted_region / all_signal}\t{frags_plotted_region}\n")

# Quantification on 4C using norm-smoothed
# First store the regions in a dictionary
quantif_regions = {}
with open(bed_file_to_quantify, 'r') as f:
    for line in f:
        if line[:5] == 'track':
            continue
        ls = line.strip().split("\t")
        quantif_regions[ls[3]] = (int(ls[1]), int(ls[2]))
# Then perform each comparison
quantif_table = table_to_quantify
with open(quantif_table, 'r') as f:
    with open(output_file_quantif_4C, 'w') as fo:
        fo.write("Sample1\tSample2\tViewpoint\tRegion\tSum_sample1\tSum_sample2\tRatio\n")
        for i, line in enumerate(f):
            if i == 0:
                # There is a header
                continue
            ls = line.strip().split("\t")
            sample1 = ls[0]
            sample2 = ls[1]
            viewpoints = ls[2].split(',')
            regions = ls[3].split(',')
            for vp in viewpoints:
                file1 = [ff for ff in os.listdir(path_4C_geo) if ff.startswith(f"{sample1}_{vp}")]
                if len(file1) != 1:
                    print(f"{sample1}_{vp}")
                    print(file1)
                    raise Exception("More than one file")
                file1 = file1[0]
                file2 = [ff for ff in os.listdir(path_4C_geo) if ff.startswith(f"{sample2}_{vp}")]
                if len(file2) != 1:
                    print(f"{sample2}_{vp}")
                    print(file2)
                    raise Exception("More than one file")
                file2 = file2[0]
                files = {}
                files[sample1] = os.path.join(path_4C_geo, file1)
                files[sample2] = os.path.join(path_4C_geo, file2)
                quantifs = {}
                # For each file
                for sample in files:
                    file = files[sample]
                    quantifs[sample] = {}
                    for region in regions:
                        quantifs[sample][region] = 0
                    with gzip.open(file, 'rb') as bdg:
                        for line_bdg in bdg:
                            line_bdg = line_bdg.decode('ascii')
                            if line_bdg[:5] == 'track':
                                continue
                            ls_bdg = line_bdg.strip().split("\t")
                            start_bdg = int(ls_bdg[1])
                            end_bdg = int(ls_bdg[2])
                            # We only consider fragments in the plotted region
                            if end_bdg < plotted_region[0]:
                                continue
                            if start_bdg > plotted_region[1]:
                                continue
                            # Then we check for each region the overlap
                            for region in regions:
                                if end_bdg > quantif_regions[region][0] and \
                                        start_bdg < quantif_regions[region][1]:
                                    quantifs[sample][region] += float(ls_bdg[3])
                for region in regions:
                    fo.write("\t".join([sample1, sample2, vp, region, \
                        str(quantifs[sample1][region]), str(quantifs[sample2][region]),
                        str(round(quantifs[sample1][region] / quantifs[sample2][region] * 100, 2))]))
                    fo.write("\n")

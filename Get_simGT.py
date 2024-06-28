#!/usr/bin/env python3
import sys
import pysam, os
from collections import defaultdict
import numpy as np

_, bed_file, outfile, *bamfiles = sys.argv

if os.path.exists(bed_file):
    in_bed = defaultdict(list)
    with open(bed_file) as data:
        for line in data:
            chrom, start, end = line.strip().split()
            start, end = int(start), int(end)
            in_bed[chrom].append(end)

else:
    print("Bedfile missing")
    sys.exit()

if os.path.exists(outfile):
    sys.exit()

sim_dict=defaultdict(dict)
pos_dict=defaultdict(list)
#print("Chr", "Start",  "End", "Base",sep="\t") #header
for bam_file in bamfiles:
    # Open the BAM file for reading
    bam = pysam.AlignmentFile(bam_file, "rb")

    snp_count=0

    # Open the BAM file using a with statement for automatic closing
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate through the records in the BAM file
        for record in bam:
            # Extract values from the specified columns
            #print(record)
            readname = record.query_name
            readnameA,readnameB=readname.split(":")[0:2]
            simpos_start,rest=readnameB.split("-")[0:2]
            simpos_end=rest.split("_")[0]
            simchr=readnameA.split("_")[3]
            simpos_start, simpos_end = int(simpos_start), int(simpos_end)
            #print(simchr, simpos_start, simpos_end, sep="\t")
            #print(in_bed[simchr])
            if any(simpos_start <= x <= simpos_end for x in in_bed[simchr]):
                    snp_count += 1
                    snp_pos = [x for x in in_bed[simchr] if x >= simpos_start and x <= simpos_end][0]
                    read_pos=len(record.query_sequence)-(simpos_end-snp_pos+1)
                    read_GT=record.query_sequence[read_pos]
                    #print(record)
                    #print(snp_pos, read_pos, read_GT, simpos_end, sep="\t")
                    #if read_GT=="N":
                    #    continue
                    #else:
                    pos_dict[snp_pos].append(read_GT)
                    sim_dict[simchr]=[pos_dict]
                        #in_bed[chrom][start:end]
                        #snp=in_bed[simchr]
                        #print(simchr, snp_pos-1,snp_pos,read_GT,snp_count, sep="\t")
#print(sim_dict)
with open(outfile, "w") as output:
    header = "\t".join(["Chr", "Start",  "End", "Base", "\n"])
    output.write(header)
    for i in sim_dict:
        #print(i)
        for posd in sim_dict[i]:
            for pos in posd:
                gts=posd[pos]
                lineinfo = "\t".join([str(i),str(pos-1),str(pos),str(gts),"\n"])
                output.write(lineinfo)
                #output.write(i,pos-1,pos,gts, sep="\t")


#!/usr/bin/env python3
import sys,re
import pysam, os
from collections import defaultdict
import numpy as np
from collections import Counter
import gzip


vcf_file1=sys.argv[1] #sim_vcf_file
vcf_file2=sys.argv[2] #source_vcf_file
individual=sys.argv[3] #target indiv index
labels=sys.argv[4] #labels for file



in_vcf_bases = defaultdict(dict)
in_vcf_gts = defaultdict(dict)  
SNPs=0  
with open(vcf_file2) as data:
    for line in data:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            fields = line.strip().split("\t")
            ind_index = int(individual)+9
            ind_ID=fields[ind_index]
            #print(ind_ID, fields[1], sep="\t")
        else:
            field_ind = line.strip().split("\t")[ind_index]
            chrom, pos,_,ref,alt = line.strip().split("\t")[0:5]
            #print(field_ind)
            genotype = field_ind.split(":")[0]
            gt1=genotype.split("/")[0] #called genotype base A
            gt2=genotype.split("/")[1] # called genotype B
            alt_list=alt.split(",")
            alt_list.insert(0,ref)
            if gt1 !="." and gt2 !=".":
                gt_base1=alt_list[int(gt1)] #called base for genotype A
                gt_base2=alt_list[int(gt2)] #called base for genotype B
                genotype_bases=[gt_base1, gt_base2]
                SNPs+=1
                in_vcf_bases[chrom][pos] = genotype_bases
                in_vcf_gts[chrom][pos] = genotype
                #print("")


def values_in_range(dictionary, lower_bound, upper_bound):
    """
    Function to find keys within a given range in a dictionary.

    Parameters:
    - dictionary: The dictionary to search through.
    - lower_bound: The lower bound of the range.
    - upper_bound: The upper bound of the range.

    Returns:
    - A keys in dictionary in set range
    """
    
    for key, value in dictionary.items():
        #return(type(key))
        # Check if the value is an integer
        if lower_bound <= int(key) <= upper_bound:
            return(key)

def count_concordant_genotypes(vcf_file1, vcf_dict_bases, vcf_dict_gts, individual): #Checks if called genotypes match "correct" genotypes
    total_sim_SNPs=0
    concord_GTs=0
    disconcord_GTs=0
    GTs_refref=0
    GTs_refalt=0
    GTs_altalt=0
    GTs_alt2alt=0
    GTs_2alt2alt=0
    GTs_ref2alt=0
    cGTs_refref=0
    dGTs_refref=0
    cGTs_refalt=0
    dGTs_refalt=0
    cGTs_altalt=0
    dGTs_altalt=0
    cGTs_alt2alt=0
    dGTs_alt2alt=0
    cGTs_ref2alt=0
    dGTs_ref2alt=0
    cGTs_2alt2alt=0
    dGTs_2alt2alt=0
    with open(vcf_file1,'rt') as simfile:
        print("labels", "Type","Concordant", "Disconcordant","Total",sep="\t")
        for line in simfile:
            if line.startswith("#"):
                continue
            else:
                chrom, pos,_,ref,alt = line.strip().split("\t")[0:5]
                field_ind=line.strip().split("\t")[9]
                #print(field_ind)
                genotypeold = field_ind.split(":")[0]
                #DP=max(field_ind.split(":")[1])
                genotype = genotypeold.replace("|", "/")
                gt1=genotype.split("/")[0] #called genotype base A
                gt2=genotype.split("/")[1] # called genotype B
                alt_list=alt.split(",")
                alt_list.insert(0,ref)
                if gt1 !="." and gt2 !=".":
                    gt_base1=alt_list[int(gt1)].upper() #called base for genotype A
                    gt_base2=alt_list[int(gt2)].upper() #called base for genotype B
                    #print(gt_base1,gt_base2)
                    #print(alt_list) # all possible bases (ref and alternative)
                    #print(ref, alt, genotype, gt1, gt2,gt_base1,gt_base2, sep="\t")
                    genotype_bases=[gt_base1, gt_base2]
                    newposcol=line.strip().split("\t")[7]
                    newstart=newposcol.split("=")[0]
                    newpos = newposcol.split("=")[1]
                    #print(newstart, newpos, pos,type(newpos), type(pos), sep="\t")
                    if pos in vcf_dict_bases[chrom]:
                        #print(pos)
                        total_sim_SNPs+=1
                        called_gtBASES=[*set([gt_base1.upper(),gt_base2.upper()])]
                        correctGTbases=vcf_dict_bases[chrom][pos]
                        correctgt=vcf_dict_gts[chrom][pos]
                        if correctgt == "0/0":
                            GTs_refref+=1
                            if sorted(genotype_bases)==sorted(correctGTbases):
                                concord_GTs+=1
                                cGTs_refref+=1
                            else:
                                #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                                disconcord_GTs+=1
                                dGTs_refref+=1       
                    else:
                        #print(newstart)
                        if newstart == "END":
                            #print(newpos)
                            gvcfSNPs=values_in_range(vcf_dict_bases[chrom], int(pos), int(newpos))
                            if hasattr(gvcfSNPs, '__len__'):
                                total_sim_SNPs+=1
                                called_gtBASES=[*set([gt_base1.upper(),gt_base2.upper()])]
                                correctGTbases=vcf_dict_bases[chrom][gvcfSNPs]
                                correctgt=vcf_dict_gts[chrom][gvcfSNPs]
                                pos=gvcfSNPs
                                if correctgt == "0/0":
                                    GTs_refref+=1
                                    concord_GTs+=1
                                    cGTs_refref+=1
                            else:
                                continue
                        else:
                            #print(pos)
                            continue
                    if correctgt == "0/1" or correctgt =="1/0":
                        GTs_refalt+=1
                        if sorted(genotype_bases)==sorted(correctGTbases):
                            concord_GTs+=1
                            cGTs_refalt+=1
                        else:
                            #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                            disconcord_GTs+=1
                            dGTs_refalt+=1
                    elif correctgt == "1/1":
                        GTs_altalt+=1
                        if sorted(genotype_bases)==sorted(correctGTbases):
                            concord_GTs+=1
                            cGTs_altalt+=1
                        else:
                            #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                            disconcord_GTs+=1
                            dGTs_altalt+=1
                    elif correctgt == "1/2" or correctgt =="2/1":
                        GTs_alt2alt+=1
                        if sorted(genotype_bases)==sorted(correctGTbases):
                            concord_GTs+=1
                            cGTs_alt2alt+=1
                        else:
                            #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                            disconcord_GTs+=1
                            dGTs_alt2alt+=1
                    elif correctgt == "0/2" or correctgt =="2/0":
                        GTs_ref2alt+=1
                        if sorted(genotype_bases)==sorted(correctGTbases):
                            concord_GTs+=1
                            cGTs_ref2alt+=1
                        else:
                            #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                            disconcord_GTs+=1
                            dGTs_ref2alt+=1
                    elif correctgt == "2/2":
                        GTs_2alt2alt+=1
                        if sorted(genotype_bases)==sorted(correctGTbases):
                            concord_GTs+=1
                            cGTs_2alt2alt+=1
                        else:
                            #print(sorted(genotype_bases),sorted(correctGTbases), correctgt, genotype,pos)
                            disconcord_GTs+=1
                            dGTs_2alt2alt+=1
                

        print(labels, "Total",concord_GTs, disconcord_GTs,total_sim_SNPs,sep="\t")
        print(labels, "0/0",cGTs_refref, dGTs_refref,GTs_refref,sep="\t")
        print(labels, "0/1",cGTs_refalt, dGTs_refalt,GTs_refalt,sep="\t")
        print(labels, "1/1",cGTs_altalt, dGTs_altalt,GTs_altalt,sep="\t")
        print(labels, "1/2",cGTs_alt2alt, dGTs_alt2alt,GTs_alt2alt,sep="\t")
        print(labels, "0/2",cGTs_ref2alt, dGTs_ref2alt,GTs_ref2alt,sep="\t")
        print(labels, "2/2",cGTs_2alt2alt, dGTs_2alt2alt,GTs_2alt2alt,sep="\t")

                    



result = count_concordant_genotypes(vcf_file1,in_vcf_bases,in_vcf_gts, individual)

# Workflow for testing genotyping and imputation methods with simulated sequencing data

This is the GWF workflow used in Zavala, Rohlfs, Moorjani, "Benchmarking for genotyping and imputation using degraded DNA for forensic applications across diverse populations",FISG 2024 (<https://doi.org/10.1016/j.fsigen.2024.103177>). 

For information on how to set up GWF please see <https://gwf.app/>

## Contact
<https://elenazavala.owlstown.net/>

---

## Parameters

The workflow described below that is used in file workflow.py iterates creates simulated data for different individuals, SNP panels, chromosomes, coverages, and DNA qualities. 

## Workflow Steps

NOTE: There are paths and options in the workflow that will need to be updated for your file system. These are noted within the workflow.py file. Within workflow.py the required software is listed for each step

1) Simulation of sequencing data with NGSNGS <https://github.com/RAHenriksen/NGSNGS>

This requires reference fasta files of the human reference genome hg38 build and vcf files with genotype data subset to the SNPs of interest. 

The output from this step is simulated sequences of three different qualities (Degraded A, Degraded B, and Modern)

2) Mapping of simulated sequences to the human reference genome (hg38) with ancient parameters

This step requires bwa <https://bio-bwa.sourceforge.net/bwa.shtml> and SAMTools <https://github.com/samtools/htslib>. 

3) Filtering of mapped bam files

Reads less than 35bps and with mapping qualities less than 25 are removed with SAMTools. 

4) Trim deamination of filtered bam files

Custom script trim_Xbps_BAM.py used to remove potential deamination from terminal 8bps by replacing these substitutions with N's

5) Genotyping of trimmed, filtered and mapped bam files

   This was performed using either ATLAS (<https://bitbucket.org/wegmannlab/atlas/wiki/Home>), GATK (<https://gatk.broadinstitute.org/hc/en-us>), or SAMTools

6) Refine genotypers from step 5

Genotype refinement was performed with either Beagle 4.1 (<https://faculty.washington.edu/browning/beagle/b4_1.html>), GATK, or BCFTtools. For Beagle 4.1 and GATK reference panels from the high coverage 1000 Genomes database were used (<https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/>). 

7) Imputation

Imputation was performed with either Beagle 5.4 (<https://faculty.washington.edu/browning/beagle/beagle.html>) or GLIMPSE2 (<https://odelaneau.github.io/GLIMPSE/>) using the same reference panels that were used for genotype refinement. 

8) Check Accuracy of genotyped SNPs

Performed by comparing genotyped SNPs in output VCF for each respective step to the original vcf used for simulating genotypes 






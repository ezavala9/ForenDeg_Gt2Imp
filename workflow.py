from gwf import Workflow, AnonymousTarget
import os
import glob
from collections import defaultdict
#gwf config set backend slurm
gwf = Workflow()


def sim_seqData(cov, panel_label, ind,chr, outputfolder):
    """Simulating single read sequence data with NGSNGS."""
    inputs = []
    outputs = [f"{outputfolder}/{panel_label}_{ind}_{chr}_DegMid2Sim_cov{cov}.fq.txt"]+ [f"{outputfolder}/{panel_label}_{ind}_{chr}_ModSim_cov{cov}.fq.txt"]+[f"{outputfolder}/{panel_label}_{ind}_{chr}_DegSim_cov{cov}.fq.txt"]#+[f"Simulated_fastq/{panel_label}_{ind}_{chr}_ModSim_cov{cov}.fq"]+[f"Simulated_fastq/{panel_label}_{ind}_{chr}_DegSim_cov{cov}.fq"]+[f"{outputfolder}/{panel_label}_{ind}_{chr}_DegMid2Sim_cov{cov}.fq"]
    options = {
        'cores': 1,
        'memory': '10g', 
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""

    module load gcc
    module load samtools

    NGSNGS/ngsngs -i hg38/fasta/hg38_chr{chr}.fa -vcf SNP_panels/{panel_label}_Chr{chr}_HG38.noindels.vcf -s 22688 -id {ind} -c {cov} -f fq -l 150 -seq se -q1 NGSNGS/Test_Examples/AccFreqL150R1.txt -o {outputfolder}/{panel_label}_{ind}_{chr}_ModSim_cov{cov}
    NGSNGS/ngsngs -i hg38/fasta/hg38_chr{chr}.fa -vcf SNP_panels/{panel_label}_Chr{chr}_HG38.noindels.vcf -s 22688 -id {ind} -c {cov} -f fq -ld norm,40,10 -seq se -q1 NGSNGS/Test_Examples/AccFreqL150R1.txt -m b7,0.024,0.36,0.68,0.0097 -o {outputfolder}/{panel_label}_{ind}_{chr}_DegSim_cov{cov}
    NGSNGS/ngsngs -i hg38/fasta/hg38_chr{chr}.fa -vcf SNP_panels/{panel_label}_Chr{chr}_HG38.noindels.vcf -s 22688 -id {ind} -c {cov} -f fq -ld norm,80,20 -seq se -q1 NGSNGS/Test_Examples/AccFreqL150R1.txt -m b7,0.01,0.52,0.42,0.0097 -o {outputfolder}/{panel_label}_{ind}_{chr}_DegMid2Sim_cov{cov}
    
    echo done>{outputfolder}/{panel_label}_{ind}_{chr}_ModSim_cov{cov}.fq.txt
    echo done > {outputfolder}/{panel_label}_{ind}_{chr}_DegSim_cov{cov}.fq.txt
    echo done > {outputfolder}/{panel_label}_{ind}_{chr}_DegMid2Sim_cov{cov}.fq.txt

    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
  

def bwa_ANCmap(input_delcheck,ref_genome, r1, bamfile, fastrim, outputfolder):
    """Mapping reads to a reference genome with `bwa` and `samtools` with ancient paraments"""
    inputs = [f"{input_delcheck}"] #used to link to sim_seqData
    outputs = [f"{outputfolder}/ANC{bamfile}_{fastrim}.bam.txt"]
    options = {
        'cores': 16,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f'''
    module load bwa
    module load samtools

    bwa aln {ref_genome} -t 16 -n 0.01 -o 2 -l 16500 {r1} > temp/ANC{bamfile}_{fastrim}.sai 
    bwa samse {ref_genome} temp/ANC{bamfile}_{fastrim}.sai {r1} >temp/ANC{bamfile}_{fastrim}.sam
    samtools view -h -b temp/ANC{bamfile}_{fastrim}.sam | samtools sort -T temp/ANC{bamfile}_{fastrim} -o {outputfolder}/ANC{bamfile}_{fastrim}.bam
    samtools index {outputfolder}/ANC{bamfile}_{fastrim}.bam
    
    rm temp/ANC{bamfile}_{fastrim}.sai
    rm temp/ANC{bamfile}_{fastrim}.sam
    rm {r1} 
    echo done> {outputfolder}/ANC{bamfile}_{fastrim}.bam.txt
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def filter_bam( input_bam_anc, outputbase,Trim, fastrim, del_base, outputfolder):
    """filter based on quality and rename RG"""
    inputs = [f"{input_bam_anc}.txt"] #used to link to bwa_ANCmap
    outputs = [f"{outputfolder}/{outputbase}_{Trim}_L35MQ25.bam.bai"]+[f"{outputfolder}/{outputbase}_{Trim}_L35MQ25.bam.txt"]
    options = {
        'cores': 1,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load samtools

    ##filtering based on length and mapping quality
    samtools view -h -q 25 {input_bam_anc}| awk 'length($10) > 35 || $1 ~ /^@/' | samtools view -bS - > temp/{outputbase}_{fastrim}_L35MQ25.bam
    
    ## Update RG for downstream processing
    samtools addreplacerg -r "@RG\tID:S1" -o {outputfolder}/{outputbase}_{Trim}_L35MQ25.bam temp/{outputbase}_{fastrim}_L35MQ25.bam
    samtools index {outputfolder}/{outputbase}_{Trim}_L35MQ25.bam
    
    ## Remove fastq and bam files no longer needed to save space
    rm -f Simulated_fastq/{del_base}.fq
    rm -f {input_bam_anc} {outputbase}_{fastrim}_L35MQ25.bam
    rm -f temp/{del_base}_{Trim}*sam temp/{del_base}_{Trim}*sai temp/{del_base}_{Trim}*bam 
    
    echo done > {outputfolder}/{outputbase}_{Trim}_L35MQ25.bam.txt
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def trimbam(inputbam, outputbam):
    """Replace putative deamination with 'N's"""
    inputs = [f"{inputbam}.txt"] #used to link to filter_bam
    outputs = [f"{outputbam}.bai"]+[f"{outputbam}.txt"]
    options = {
        'cores': 1,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f'''
    module load samtools 

    python trim_Xbps_BAM.py {inputbam} {outputbam}
    samtools index {outputbam}

    echo done > {outputbam}.txt
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

  
def atlas_geno_bam(input_bam, input_base, pan,chr):
    """genotype bams with ATLAS"""
    inputs = [f"{input_bam}.txt"] #used to link to trim bam
    outputs = [f"GT_ATLAS/{input_base}_MaximumLikelihood_subset{pan}.vcf"]
    options = {
        'cores': 1,
        'memory': '70g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load gcc
    module load bcftools
    module load samtools
    module load openblas
    module load lapack
    module load bedtools
    
    ## Genotyping with ATLAS
    /global/scratch/p2p3/pl1_moorjani/ezavala9/dev/atlas/atlas task=call method=MLE infoFields=DP formatFields=GT,AD,DP,PL bam={input_bam}  fasta=/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg38/fasta/hg38.fa  out=GT_ATLAS/At_GT_{input_base}
    
    ## Subset to SNP panel of interest
    bedtools intersect -header -wa -a GT_ATLAS/At_GT_{input_base}_MaximumLikelihood.vcf.gz -b SNP_panels/{pan}_chr{chr}_noindels.bed > GT_ATLAS/{input_base}_MaximumLikelihood_subset{pan}.vcf  
    
    #Remove previous genotype file to save space
    rm GT_ATLAS/At_GT_{input_base}_MaximumLikelihood.vcf.gz 
    
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_geno_bam(input_bam, input_base, pan,ref,chr):
    """genotype bams SAMTools"""
    inputs = [f"{input_bam}.txt"]
    outputs = [f"GT_Samtools/{input_base}_Samtools_subset{pan}.vcf" ]
    options = {
        'cores': 1,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load bcftools
    module load samtools
    module load bedtools

    ##Genotype with SAMTools
    samtools mpileup -C50 --redo-BAQ --min-BQ 25 --output-tags DP,AD -f {ref} --BCF {input_bam} | bcftools call --consensus-caller -Ov > GT_Samtools/{input_base}_all.vcf
    
    ## Subset to SNP panel of interest
    bedtools intersect -header -wa -a GT_Samtools/{input_base}_all.vcf -b SNP_panels/{pan}_chr{chr}_noindels.bed > GT_Samtools/{input_base}_Samtools_subset{pan}.vcf  
    
    #Remove previous genotype file to save space
    rm GT_Samtools/{input_base}_all.vcf 
    
   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
 
def GATK_geno_bam(input_bam, input_base, pan,ref,chr):
    """genotype bams with GATK"""
    inputs = [f"{input_bam}.txt"]
    outputs = [f"GT_GATK/{input_base}_GATK_subset{pan}.vcf" ]
    options = {
        'cores': 4,
        'memory': '15g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load java/1.8.0_121
    module load samtools
    module load bedtools
    module load bcftools

    ## Reheader bamfile for input into GATK
    samtools view -H {input_bam} | sed 's,^@RG.*,@RG\tID:S1\tSM:S1\tLB:None\tPL:Illumina,g' | samtools reheader - {input_bam} > temp/{input_base}_GATK_Hfix.bam
    samtools index temp/{input_base}_GATK_Hfix.bam

    ##Genotype with GATK
    /global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/genome_analysis/gatk-4.1.9.0/gatk --java-options "-Xmx10g" HaplotypeCaller -R {ref} -I temp/{input_base}_GATK_Hfix.bam -O GT_GATK/{input_base}.g.vcf.gz -ERC GVCF
    
    ## Subset to SNP panel of interest
    bedtools intersect -header -wa -a GT_GATK/{input_base}.g.vcf.gz -b /global/scratch/p2p3/pl1_moorjani/ezavala9/forensic/2023_NGSNGS_sims/SNP_panels/{pan}_chr{chr}_noindels.bed | bcftools view  -i  'MIN(FMT/DP)>0'> GT_GATK/{input_base}_GATK_subset{pan}.vcf  
    
    #Remove previous genotype file to save space
    rm -f temp/{input_base}_GATK_Hfix.bam  
    
   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GTrefine_BG(input_vcf, outfile,reffile, pan,chr):
    """Refine genotypes with Beagle 4.1"""
    inputs = [{input_vcf}] # used to link to output of ATLAS
    outputs = [f"Refine_Bg/{outfile}.vcf.gz" ]+[f"Refine_Bg/{outfile}.gp0.99.vcf.gz"]
    options = {
        'cores': 4,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load java
    module load bcftools

    ## Refine with Beagle 4.1 
    java -Xss10m -Xmx16g -jar beagle.27Jan18.7e1.jar gl={input_vcf} ref=Ref_1000G/{reffile}_450subset_chr{chr}.vcf.gz out=Refine_Bg/{outfile} impute=false nthreads=4 map=plink.GRCh38.map/plink.chr{chr}.GRCh38.map gprobs=true 

    ##Subset to SNPs with a probabilit of at least 99%
	bcftools view -i 'MAX(GP[*])>=0.99' Refine_Bg/{outfile}.vcf.gz  -Oz -o Refine_Bg/{outfile}.gp0.99.vcf.gz 
    
   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GTrefine_GATK(input_vcf, outfile,reffile, pan,chr):
    """Refine genotypes with GATK"""
    inputs = [{input_vcf}] #used to link to output of GATK genotyping
    outputs = [f"Refine_GATK/{outfile}.gp0.99.vcf.gz" ]+[f"Refine_GATK/{outfile}.vcf.gz" ]
    options = {
        'cores': 4,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load java/1.8.0_121
    module load bcftools

    ## Refine with GATK
    /global/scratch/p2p3/pl1_moorjani/SHARED_LAB/bin/genome_analysis/gatk-4.1.9.0/gatk --java-options "-Xmx4g" CalculateGenotypePosteriors -V {input_vcf} -O Refine_GATK/{outfile}.vcf.gz -supporting Ref_1000G/{reffile}_450subset_chr{chr}.vcf.gz 

    ##Subset to SNPs with a probabilit of at least 99%
	bcftools view -i'FMT/GQ>20' Refine_GATK/{outfile}.vcf.gz -Oz -o Refine_GATK/{outfile}.gp0.99.vcf.gz 
    
   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GTrefine_ST(input_vcf, outfile):
    """Refine genotypes with bcftools"""
    inputs = [{input_vcf}]
    outputs = [f"Refine_ST/{outfile}.gp0.99.vcf.gz" ]
    options = {
        'cores': 4,
        'memory': '10g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load bcftools

    ##Subset to SNPs with a probabilit of at least 99%
    bcftools view -i '%QUAL>=20' {input_vcf} -Oz -o Refine_ST/{outfile}.gp0.99.vcf.gz 
    
   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf2vcf_check(code, outbase, sim_vcf, ref_vcf, ind,labels,outfolder):
    """genotype bams with different references"""
    inputs = [sim_vcf]
    outputs = [f"{outfolder}/ANC_{outbase}.txt"]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    python {code} {sim_vcf} {ref_vcf} {ind} {labels} > {outfolder}/ANC_{outbase}.txt
   
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def Impute_BG(input_vcf, outfile,reffile, pan,chr):
    """Impute with Beagle 5.4"""
    inputs = [{input_vcf}] #links to output to of beagle 4.1 genotype refinement
    outputs = [f"Impute_Bg/{outfile}.gp0.99_subset{pan}_IMPonly.txt" ]+[f"Impute_Bg/{outfile}.gp0.99_subset{pan}.vcf" ]
    options = {
        'cores': 4,
        'memory': '25g',
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load java
    module load bcftools
    module load bedtools

    ## Imputation with Beagle 5.4
    java -Xss10m -Xmx16g -jar /global/scratch/p2p3/pl1_moorjani/ezavala9/dev/beagle.22Jul22.46e.jar gt={input_vcf} ref=Ref_1000G/{reffile}_450subset_chr{chr}.vcf.gz out=Impute_Bg/{outfile} impute=true nthreads=4 map=/global/scratch/p2p3/pl1_moorjani/ezavala9/dev/plink.GRCh38.map/plink.chr{chr}.GRCh38.map gp=true 
    
    ## Subset to genotypes with at least 99% genotyper probability
	bcftools view -i 'MAX(GP[*])>=0.99' Impute_Bg/{outfile}.vcf.gz  -Oz -o Impute_Bg/{outfile}.gp0.99.vcf.gz 

    ## Subset to SNPs of interest
    bedtools intersect -header -wa -a Impute_Bg/{outfile}.gp0.99.vcf.gz -b SNP_panels/{pan}_chr{chr}_noindels.bed > Impute_Bg/{outfile}.gp0.99_subset{pan}.vcf
    
    ## Subset to only imputed SNPs
    grep IMP Impute_Bg/{outfile}.gp0.99_subset{pan}.vcf > Impute_Bg/{outfile}.gp0.99_subset{pan}_IMPonly.txt

   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def Impute_Glimpse(input_bam, outfile,ref, pan,chr):
    """Impute genotypes with Glimpse2 from Samtools data"""
    inputs = [input_bam] # link to refined genotype outputs
    outputs = [f"Impute_Glimpse4Mb_ST/{outfile}_ligated.bcf" ] 
    options = {
        'cores': 4,
        'memory': '25g', 
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load gcc
    module load bedtools
    module load bcftools

    bcftools index -f {input_bam}
    #bcftools +tag2tag {input_bam} -- -r --gp-to-gl | bcftools view -Oz -o {input_bam}_GL.vcf.gz
    #bcftools +tag2tag {input_bam}_GL.vcf.gz -- -r --gl-to-pl | bcftools view -Oz -o {input_bam}_PL.vcf.gz

    #REF=Ref_GLIMPSE/Bins4Mb/{ref}450_chr{chr}.split
    #BAM={input_bam}

    #bcftools index -f {input_bam}_PL.vcf.gz

    while IFS="" read -r LINE || [ -n "$LINE" ]; 
    do   
        printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)
        CHR=$(echo $LINE | cut -d" " -f2)
        REGS=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f1)
        REGE=$(echo $IRG | cut -d":" -f 2 | cut -d"-" -f2)
        OUT=Impute_Glimpse4Mb_ST/Chunks/{outfile}
        
        glimpse/phase/bin/GLIMPSE2_phase --input-gl {input_bam} --reference Ref_GLIMPSE/Bins4Mb/{ref}450_chr{chr}.split_"$CHR"_"$REGS"_"$REGE".bin --output "$OUT"_"$CHR"_"$REGS"_"$REGE".bcf
    done < Ref_GLIMPSE/Bins4Mb/{ref}chunks.chr{chr}.txt

    #ligate piece together
    LST=Impute_Glimpse4Mb_ST/Chunks/list.{outfile}.txt
    ls -1v Impute_Glimpse4Mb_ST/Chunks/{outfile}_*.bcf > $LST

    OUT=Impute_Glimpse4Mb_ST/{outfile}_ligated.bcf
    glimpse/ligate/bin/GLIMPSE2_ligate --input $LST --output $OUT

   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def Impute_Glimpse_subset(outfile,ref, pan,chr,cutoff):
    """Subset Imputed genotypes with Glimpse2 for output from GLIMPSE with Beagle 4.1 output"""
    inputs = [f"Impute_Glimpse4Mb/{outfile}_ligated.bcf" ]
    outputs = [f"Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.gp{cutoff}.vcf.gz" ] +[f"Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.vcf.gz" ] +[f"Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.gp0.99_imponly.vcf.gz"]
    options = {
        'cores': 4,
        'memory': '25g', 
        'walltime': '72:00:00',
        'queue': '', ### UPDATE DEPENDENT ON SYSTEM
        'account': '', ### UPDATE DEPENDENT ON SYSTEM
    }

    spec = f"""
    module load gcc
    module load bedtools
    module load bcftools

    #convert bcf to vcf.gz
    bcftools convert -Oz -o Impute_Glimpse4Mb/{outfile}_ligated.vcf.gz Impute_Glimpse4Mb/{outfile}_ligated.bcf
    
    #subset to panel
    bedtools intersect -header -wa -a Impute_Glimpse4Mb/{outfile}_ligated.vcf.gz -b SNP_panels/{pan}_chr{chr}_noindels.bed | bcftools view -Oz -o Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.vcf.gz
    
    #subsetp to gp0.99
    bcftools view -i 'MAX(GP[*])>={cutoff}' -Ov -o Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.gp{cutoff}.vcf.gz Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.vcf.gz 

    #subset to imp only 
    bcftools view -i 'DS[*]!=0' Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.gp0.99.vcf.gz | bcftools view -i 'DS[*]!=1' | bcftools view -i 'DS[*]!=2' -Ov -o Impute_Glimpse4Mb/{outfile}_ligated_subset{pan}only.gp0.99_imponly.vcf.gz


   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def Impute_GlimpseST_subset(outfile,ref, pan,chr,cutoff):
    """Subset Imputed genotypes with Glimpse2 for output from GLIMPSE with SAMTOOLS refinement"""
    inputs = [f"Impute_Glimpse4Mb_ST/{outfile}_ligated.bcf" ]
    outputs = [f"Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.gp{cutoff}.vcf.gz" ] +[f"Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.vcf.gz" ] +[f"Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.gp0.99_imponly.vcf.gz"]
    options = {
        'cores': 4,
        'memory': '25g', 
        'walltime': '72:00:00',
        'queue': 'savio3_htc',
        'account': 'co_moorjani',
    }

    spec = f"""
    module load gcc
    module load bedtools
    module load bcftools

    #convert bcf to vcf.gz
    bcftools convert -Oz -o Impute_Glimpse4Mb_ST/{outfile}_ligated.vcf.gz Impute_Glimpse4Mb_ST/{outfile}_ligated.bcf
    
    #subset to panel
    bedtools intersect -header -wa -a Impute_Glimpse4Mb_ST/{outfile}_ligated.vcf.gz -b SNP_panels/{pan}_chr{chr}_noindels.bed | bcftools view -Oz -o Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.vcf.gz
    
    #subsetp to gp0.99
    bcftools view -i 'MAX(GP[*])>={cutoff}' -Ov -o Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.gp{cutoff}.vcf.gz Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.vcf.gz 

    #subset to imp only 
    bcftools view -i 'DS[*]!=0' Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.gp0.99.vcf.gz | bcftools view -i 'DS[*]!=1' | bcftools view -i 'DS[*]!=2' -Ov -o Impute_Glimpse4Mb_ST/{outfile}_ligated_subset{pan}only.gp0.99_imponly.vcf.gz


   """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



### column indexs for the individuals in the genotype VCF used for simulations
base_10indiv=[2309,2312,2866, 2878, 2855, 2942, 2519, 2524, 3146, 3181]
### Panel IDs 
panel=["FORCE", "TRI" , "25K", "95K", "HO"]


for ind in base_10indiv:
    for  pan in panel: 
        for chr in [1,14]: 
            for cov in [0.1,0.5,1,5,10]: 
##                simulate and map reads for various coverages
                gwf.target_from_template(name=f"sim{pan}files_{ind}_{chr}_{cov}",template=sim_seqData(cov=cov,panel_label=pan, ind=ind, chr=chr, outputfolder=f"Simulated_fastq"))
                for qual in ["Mod", "Deg", "DegMid2"]:                    
                    ##map fastqs 
                    gwf.target_from_template(name=f"ANCmap_{qual}_{pan}_{ind}_{chr}cov{cov}",template=bwa_ANCmap(ref_genome=name, input_delcheck=f'Simulated_fastq/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}.fq.txt',r1=f'Simulated_fastq/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}.fq', bamfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}',fastrim=f"NoTrim", outputfolder=f"Simulated_mapped_bams"))
                   # # Filter bams for genotyping (L35MQ25)
                    gwf.target_from_template(name=f"filter_notrim_{qual}_{pan}_{ind}_{chr}cov{cov}",template=filter_bam(input_bam_anc=f'Simulated_mapped_bams/ANC{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_NoTrim.bam',  outputbase=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}', Trim=f"NoTrim", fastrim=f"NoTrim", del_base=f"{pan}_{ind}_{chr}_{qual}Sim_cov{cov}", outputfolder=f"FilterL35MQ25"))               
                 #   Trim mapped filtered untrimmed Bams
                    gwf.target_from_template(name=f"ANCTrimMapFilt_{qual}_{pan}_{ind}_{chr}cov{cov}",template=trimbam(inputbam=f'FilterL35MQ25/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_NoTrim_L35MQ25.bam', outputbam=f'FilterL35MQ25/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_BamTrim_L35MQ25.bam'))
                    for Trim in ['BamTrim']:
                    # ### Atlas Genotyping and concord check
                        gwf.target_from_template(name=f"ATgeno_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=atlas_geno_bam(input_bam=f'FilterL35MQ25/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_L35MQ25.bam',input_base=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}', pan=pan, chr=chr))
                        gwf.target_from_template(name=f"GTvcf2vcfConcord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}", sim_vcf=f"GT_ATLAS/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_MaximumLikelihood_subset{pan}.vcf", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}", outfolder=f"Check_VCF_GT"))
                        for Ref in ['Afr', 'Amer', 'EA', 'Kg']:
                            gwf.target_from_template(name=f"GTrefine_BG{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=GTrefine_BG(input_vcf=f'GT_ATLAS/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_MaximumLikelihood_subset{pan}.vcf', outfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',reffile=Ref, chr=chr,pan=pan))
                            gwf.target_from_template(name=f"GTvcf2vcfRefineConcord99_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=vcf2vcf_check(code=f"Vcf2vcfgz_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99", sim_vcf=f"Refine_Bg/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}.gp0.99.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99", outfolder=f"Check_VCF_refine"))
                    ## Samtools Genotyping 
                        gwf.target_from_template(name=f"STgeno_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=samtools_geno_bam(input_bam=f'FilterL35MQ25/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_L35MQ25.bam',input_base=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}', pan=pan, ref=f"hg38/fasta/hg38.fa", chr=chr))
                        gwf.target_from_template(name=f"GTvcf2vcfConcord_ANC_ST_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}", sim_vcf=f"GT_Samtools/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_Samtools_subset{pan}.vcf", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_ST_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}",outfolder=f"Check_VCF_GT")) 
                    ## Samtools Refine
                        gwf.target_from_template(name=f"GTrefine_ST_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=GTrefine_ST(outfile=f"ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}", input_vcf=f"GT_Samtools/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_Samtools_subset{pan}.vcf"))
                        gwf.target_from_template(name=f"GTvcf2vcfRefineConcord_ANC_ST_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=vcf2vcf_check(code=f"Vcf2vcfgz_concordSummary.py", outbase=f"ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_gp99", sim_vcf=f"Refine_ST/ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}.gp0.99.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_ST_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_gp99",outfolder=f"Check_VCF_refine"))
                     # GATK Genotyping 
                        gwf.target_from_template(name=f"GATKgeno_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=GATK_geno_bam(input_bam=f'FilterL35MQ25/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_L35MQ25.bam',input_base=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}', pan=pan, ref=f"/global/scratch/p2p3/pl1_moorjani/SHARED_LAB/DATASETS/resources/hg38/fasta/hg38.fa", chr=chr))
                        gwf.target_from_template(name=f"GTvcf2vcfConcord_ANC_GATK_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"GATK_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}", sim_vcf=f"GT_GATK/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_GATK_subset{pan}.vcf", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_GATK_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}",outfolder=f"Check_VCF_GT"))
                    ## GATK refinement
                        for Ref in ['Afr', 'Amer', 'EA', 'Kg']:
                            gwf.target_from_template(name=f"GTrefine_GATK{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=GTrefine_GATK(input_vcf=f'GT_GATK/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_GATK_subset{pan}.vcf', outfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_GATK',reffile=Ref, chr=chr,pan=pan))
                            gwf.target_from_template(name=f"GTvcf2vcfRefineConcord99_GATK_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=vcf2vcf_check(code=f"Vcf2vcfgz_concordSummary.py", outbase=f"GATK_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99", sim_vcf=f"Refine_GATK/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_GATK.gp0.99.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_GATK_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99",outfolder=f"Check_VCF_refine"))

 

###Imputation for only larger panels with atleast 0.5 coverage

    for  pan in ["25K", "95K", "HO"]:
        for chr in [1,14]:   
            for cov in [0.5,1,5,10]: #  
                for qual in ["Mod", "Deg", "DegMid2"]:  
                    for Trim in ['BamTrim']:
                        for Ref in ['Afr', 'Amer', 'EA', 'Kg']:    
                            ##Imputation with Beagle
                            gwf.target_from_template(name=f"GTIMPUTE_BG{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=Impute_BG(input_vcf=f'Refine_Bg/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}.gp0.99.vcf.gz', outfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',reffile=Ref, chr=chr,pan=pan))
                            gwf.target_from_template(name=f"GTvcf2vcfIMPUTEConcord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gp99all", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99_IMPall", sim_vcf=f"Impute_Bg/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}.gp0.99_subset{pan}.vcf", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99_IMPall", outfolder=f"Check_VCF_IMP"))
                            gwf.target_from_template(name=f"GTvcf2vcfIMPUTEConcord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gp99ImpOnly", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99_IMP", sim_vcf=f"Impute_Bg/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}.gp0.99_subset{pan}_IMPonly.txt", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99_IMP", outfolder=f"Check_VCF_IMP"))
                            # ##Impute with GLIMPSE2 current for 4Mb
                            gwf.target_from_template(name=f"GTIMPUTE_Glimpse2{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=Impute_Glimpse(input_bam=f'Refine_Bg/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}.gp0.99.vcf.gz', outfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',ref=Ref, chr=chr,pan=pan))
                            gwf.target_from_template(name=f"GTIMPUTEsubset_Glimpse2{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=Impute_Glimpse_subset(outfile=f'{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',ref=Ref, chr=chr,pan=pan,cutoff=0.99))
                            gwf.target_from_template(name=f"GTvcf2vcfIMPUTEGlimpse_Concord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gp99", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99_IMPglimpse4Mb", sim_vcf=f"Impute_Glimpse4Mb/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_ligated_subset{pan}only.gp0.99.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99_IMPglimpse", outfolder=f"Check_VCF_IMP_4Mb"))
                            gwf.target_from_template(name=f"GTvcf2vcfIMPUTEGlimpse_Concord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gpIMP", template=vcf2vcf_check(code=f"Vcf2vcfgz_concordSummary.py", outbase=f"AT_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_IMPglimpseAll4Mb", sim_vcf=f"Impute_Glimpse4Mb/{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_ligated_subset{pan}only.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_AT_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_IMPglimpseIMP", outfolder=f"Check_VCF_IMP_4Mb"))

                            ###SAMTOOLS with Glimpse
                            gwf.target_from_template(name=f"ST_GTIMPUTE_Glimpse2{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=Impute_Glimpse(input_bam=f'Refine_ST/ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}.gp0.99.vcf.gz', outfile=f'ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',ref=Ref, chr=chr,pan=pan))
                            gwf.target_from_template(name=f"ST_GTIMPUTEsubset_Glimpse2{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}", template=Impute_GlimpseST_subset(outfile=f'ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}',ref=Ref, chr=chr,pan=pan,cutoff=0.99))
                            gwf.target_from_template(name=f"ST_GTvcf2vcfIMPUTEGlimpse_Concord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gp99", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_gp99_IMPglimpse4Mb", sim_vcf=f"Impute_Glimpse4Mb_ST/ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_ligated_subset{pan}only.gp0.99.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_ST_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_gp99_IMPglimpse", outfolder=f"Check_VCF_IMP_4Mb_ST"))
                            gwf.target_from_template(name=f"ST_GTvcf2vcfIMPUTEGlimpse_Concord_ANC_{qual}_{pan}_{ind}_{chr}cov{cov}{Trim}_{Ref}gpIMPonly", template=vcf2vcf_check(code=f"Vcf2vcf_concordSummary.py", outbase=f"ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_IMPglimpseIMP4Mb", sim_vcf=f"Impute_Glimpse4Mb_ST/ST_{pan}_{ind}_{chr}_{qual}Sim_cov{cov}_{Trim}_{Ref}_ligated_subset{pan}only.gp0.99_imponly.vcf.gz", ref_vcf=f"SNP_panels/{pan}_Chr{chr}_HG38.noindels.vcf", ind=ind,labels=f"ANC_ST_{pan}_{qual}_{ind}_{chr}_{cov}_{Trim}_{Ref}_IMPglimpseIMP", outfolder=f"Check_VCF_IMP_4Mb_ST"))

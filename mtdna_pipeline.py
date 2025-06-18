
#!/usr/bin/env python

import sys
modules_path = '/home/gencardio/tfm_eli/PIPELINE/modules'
sys.path.append(modules_path)

from  trimming import trimming_fastqs, fastp_metrics, eliminate_duplicates
from fastqc import fastqc_trimmed
from mapping import download_reference_genome, mapping_fastqs, from_sam_to_bam, from_unsorted_to_sorted_bam, index_bam, download_reference_chromosome, remove_files, sorted_bam_files
from metrics import  percent_duplications, Create_sequence_dictionary, obtain_interval_list,  CollectHsMetrics, collect_insert_size_metrics, collect_aligment_summary_metrics, obtain_metrics_csv
from variant_calling import GATK_module,  Get_pileup_summaries


if __name__ == "__main__":

    # trimming_fastqs

    input_dir = "/home/gencardio/tfm_eli/Mostres_30_ultimes"
    output_dir = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"

    trimming_fastq = trimming_fastqs(input_dir, output_dir)

    #fastqc trimmed
    input_dir_fastqs = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"

    fastqc = fastqc_trimmed(input_dir_fastqs)

    #fastp_metrics

    input_dir_json = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"

    fastqs_metrics = fastp_metrics(input_dir_json)

    # #Download reference genome hg38

    genome_version = 'hg38'
    download_genome_dir = f"/home/gencardio/tfm_eli/PIPELINE/{genome_version}"

    genome_hg38 = download_reference_genome(download_genome_dir,genome_version)

    #Obtain .sam files hg38

    input_fastq = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"
    genome_reference = 'hg38'
    final_sorted_bam_files_hg38 = sorted_bam_files(input_fastq, genome_reference)
    mapping_fastq_hg38 = mapping_fastqs(input_fastq, trimming_fastq, genome_hg38,genome_reference)

    #Convert files from sam to bam, obtain unsorted bam files
    unsorted_bam_files_hg38 = from_sam_to_bam(mapping_fastq_hg38,final_sorted_bam_files_hg38)

    #Obtain sorted bam files
    sorted_bam_files_hg38 = from_unsorted_to_sorted_bam(unsorted_bam_files_hg38)

    #Remove duplicates with picard
    picard_path = "/home/gencardio/tfm_eli/PIPELINE/picard.jar"
    remove_duplicates_hg38 = eliminate_duplicates(sorted_bam_files_hg38, picard_path)

    #Index sorted bam files
    index_bam_files_hg38 = index_bam(remove_duplicates_hg38)

    #Eliminate .sam files
    remove_sam_files_hg38 = remove_files(index_bam_files_hg38)

##Metrics

    #duplications %

    percent_duplications_metrics = percent_duplications(index_bam_files_hg38, fastqs_metrics)

    #Collect Insert Size metrics

    insert_size_metrics = collect_insert_size_metrics(remove_duplicates_hg38, picard_path, percent_duplications_metrics)

     #Collect Alingment Summary Metrics

    alignment_summary_metrics = collect_aligment_summary_metrics(remove_duplicates_hg38, picard_path, genome_hg38, insert_size_metrics)

   
    #CollectHsMetrics

    #1. Obtain hg38 dictionary


    sequence_dictionary = Create_sequence_dictionary(picard_path, genome_hg38)
    
    #2. bed to interval.list

    bed_file_path = '/home/gencardio/tfm_eli/chrM/chrM'
    interval_list = obtain_interval_list(picard_path, sequence_dictionary, bed_file_path)


    #3. Obtain HsMetrics

    Hs_metrics = CollectHsMetrics(remove_duplicates_hg38, picard_path, genome_hg38, interval_list, alignment_summary_metrics)

    Metrics_excel = obtain_metrics_csv(Hs_metrics, output_dir)




    
    




    # #Obtain reference files for GATK


    # #GATK 

    # gatk_version = '4.6.0.0'
    # gatk_container_path = '~/tfm_eli'
    # reference_path = '/home/gencardio/tfm_eli/hg38'
    # GATK_module_hg38 = GATK_module(remove_duplicates_hg38,gatk_container_path, gatk_version, genome_reference, reference_path)

    #Get_pileup_summaries

    # common_biallelic_path = '/home/gencardio/tfm_eli/chrM/gnomad.genomes.v3.1.sites.chrM.vcf.bgz'
    # get_pileup_summaries_hg38 = Get_pileup_summaries(remove_duplicates_hg38, common_biallelic_path,gatk_version)



    # #Download mithocondrial chromosome:
    
    # genome_version = 'hg38'
    # chromosome_name = 'chrM'
    # download_genome_dir = f"/home/gencardio/tfm_eli/PIPELINE/{genome_version}/{chromosome_name}"
    # chrM_hg38 = download_reference_chromosome(download_genome_dir,genome_version,chromosome_name)

    # #Obtain .sam files chrM
    # input_fastq = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"
    # genome_reference = 'chrM'
    # final_sorted_bam_files_chrM = sorted_bam_files(input_fastq, genome_reference)
    # mapping_fastq_chrM = mapping_fastqs(input_fastq, trimming_fastq, chrM_hg38, genome_reference)

    # #Convert files from sam to bam, obtain unsorted bam files
    # unsorted_bam_files_chrM = from_sam_to_bam(mapping_fastq_chrM,final_sorted_bam_files_chrM)

    # #Obtain sorted bam files
    # sorted_bam_files_chrM = from_unsorted_to_sorted_bam(unsorted_bam_files_chrM)
    
    # #Index sorted bam files
    # index_bam_files_chrM = index_bam(sorted_bam_files_chrM)

    # #Eliminate .sam files
    # remove_sam_files_chrM = remove_files(index_bam_files_chrM)

    # #Remove duplicates with picard
    # picard_path = "/home/gencardio/tfm_eli/PIPELINE/picard.jar"
    # remove_duplicates_chrM = eliminate_duplicates(sorted_bam_files_chrM, picard_path)

    # #Obtain reference files for GATK
    # reference_path = '/home/gencardio/tfm_eli/chrM'
    # Sequence_dictionary = Create_sequence_dictionary(picard_path, reference_path, genome_reference)

    # #GATK 
    # gatk_version = '4.6.0.0'
    # gatk_container_path = '~/tfm_eli'
    # GATK_module_chrM = GATK_module(remove_duplicates_chrM,gatk_container_path, gatk_version,genome_reference, reference_path)
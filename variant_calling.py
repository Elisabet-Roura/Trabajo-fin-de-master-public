#!/usr/bin/env python

import sys
import os
import subprocess
import logging


#Pasos que estic seguint del docker 







def GATK_module(removed_duplicates_files_list,gatk_container_path, gatk_version, genome_reference, reference_path):

    vcf_files_list = []
# Workflow Mitochondrial short variant discovery (SNVs + Indels) https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels
    for removed_duplicates_file in removed_duplicates_files_list:

        variant_calling_path = removed_duplicates_file.split(genome_reference)[0]
        variant_calling_path = variant_calling_path.replace('mapping','variant_calling') + genome_reference
        # print(variant_calling_path)

        if not os.path.exists(variant_calling_path):
            os.makedirs(variant_calling_path)
            
    
        BAM_FOLDER = variant_calling_path.replace('variant_calling', 'mapping')
        READY_BAM_NAME = os.path.basename(removed_duplicates_file)

        VCF_FOLDER = variant_calling_path
        MUTECT2_VCF_NAME = READY_BAM_NAME.replace('.sorted.removed_duplicates.bam', '.mutect2.vcf.gz')

        REF_DIRNAME = reference_path
        #REF_DIRNAME = genome_path.split(f'/{genome_reference}.fa')[0]
        REF_FASTA_NAME = f'{genome_reference}.fa'
        sample_name = READY_BAM_NAME.split('_')[0]

        vcf_final_path = os.path.join(VCF_FOLDER,MUTECT2_VCF_NAME)


    
#info = https://gatk.broadinstitute.org/hc/en-us/articles/360035889991--How-to-Run-GATK-in-a-Docker-container
#mirar versio docker --version
#Get the GATK container imagedocker pull broadinstitute/gatk:4.6.0.0
# per entrar en el container docker run -v ~/tfm_eli:/gatk/data -it broadinstitute/gatk:4.6.0.0


        if not os.path.exists(vcf_final_path):

            cmd_mutect2 = f'docker run -v {BAM_FOLDER}:/bam_data/ -v {VCF_FOLDER}:/vcf_data/ -v {REF_DIRNAME}:/bundle/ -it broadinstitute/gatk:{gatk_version} gatk Mutect2 -R /bundle/{REF_FASTA_NAME} -L chrM   --mitochondria-mode -I /bam_data/{READY_BAM_NAME} -O /vcf_data/{MUTECT2_VCF_NAME}  '
            # cmd_mutect2 = f'docker run -v "${BAM_FOLDER}:/bam_data/" -v "${VCF_FOLDER}:/vcf_data/" -v "${REF_DIRNAME}:/bundle/" -it "broadinstitute/gatk:latest gatk Mutect2"  -I "/bam_data/${READY_BAM_NAME}" -O "/vcf_data/${MUTECT2_VCF_NAME}" -L "chrM" -R "/bundle/${REF_FASTA_NAME}"'
            p1 = subprocess.run(cmd_mutect2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(cmd_mutect2)
            

            #output = p1.stdout.decode("UTF-8")
            error = p1.stderr.decode("UTF-8")
            if p1.returncode != 0:
                msg = f" ERROR: Could not run Mutectc2 for sample {sample_name}"
                logging.error(msg)
                logging.error(error)
                sys.exit()
            
            if p1.returncode == 0:
                vcf_files_list.append(vcf_final_path)
        
        if os.path.exists(vcf_final_path):
            vcf_files_list.append(vcf_final_path)
    
    print(vcf_files_list)
    
    return vcf_files_list


def Get_pileup_summaries(removed_duplicates_files_list,common_biallelic,gatk_version):
   
    pileup_table_list = []
    common_biallelic_name = os.path.basename(common_biallelic)
    indexed_common_biallelic = common_biallelic.replace('vcf.bgz','vcf.bgz.tbi')
    common_biallelic_path = common_biallelic.split(common_biallelic_name)[0]

    # print(sample_name,vcf_file_name,vcf_file_path, common_biallelic_name, common_biallelic_path) 

    print(indexed_common_biallelic)
    if not os.path.exists(indexed_common_biallelic):
        
        cmd_index = f'docker run -v {common_biallelic_path}:/common_biallelic_data/\
        -it broadinstitute/gatk:{gatk_version} gatk IndexFeatureFile\
        -I /common_biallelic_data/{common_biallelic_name}'
        subprocess.run(cmd_index, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(cmd_index)
    
    for removed_duplicates_file in removed_duplicates_files_list:

        bam_file = os.path.basename(removed_duplicates_file)
        sample_name = bam_file.split('_')[0]
        bam_file_path = removed_duplicates_file.split(bam_file)[0]
        pileup_table_path = bam_file_path.replace('mapping', 'variant_calling')
        pileup_file_name = sample_name + '.pileup.table'
        pileup_table_file = os.path.join(pileup_table_path, pileup_file_name)


        print(bam_file, sample_name, bam_file_path, pileup_table_path, pileup_table_file)
        
        if not os.path.exists(pileup_table_file):
            cmd_get_pileup_summaries=f'docker run  -v {bam_file_path}:/bam_data/\
            -v {common_biallelic_path}:/common_biallelic_data/\
            -v {pileup_table_path}:/pileups_table/ -it broadinstitute/gatk:{gatk_version} gatk GetPilupSummaries\
            -I /bam_data/{bam_file}\
            -V /common_biallelic_data/{common_biallelic_name}\
            -L common_biallelic_data/{common_biallelic_name}\
            -O /pileups_table/pileups.table'
            p1=subprocess.run(cmd_get_pileup_summaries, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            error = p1.stderr.decode("UTF-8")
            if p1.returncode != 0:
                msg = f" ERROR: Could not run Lancet for sample {sample_name}"
                logging.error(msg)
                logging.error(error)
                sys.exit()

            if p1.returncode == 0:
                pileup_table_list.append(pileup_table_file)
        
        if os.path.exists(pileup_table_file):
            pileup_table_list.append(pileup_table_file)
    
    return pileup_table_list


# def calculate_contamination ():

#     cmd_contamination = f'docker run -v
#      gatk CalculateContamination \
#    -I pileups.table \
#    -O contamination.table'


def filter_mutect_calls(vcf_files_list, genome_path, reference_path):

    

    for vcf_file in vcf_files_list: 

        vcf_file_name = os.path.basename(vcf_file)
        vfc_file_path = vcf_file.split(vcf_file_name)[0]



        reference_file = f'{genome_reference}.fa'
        sample_name = vcf_file_name.split('_')[0]
        filtered_vcf = sample_name + '_' + genome_reference +'_filtered.vcf_gz'



#         if not os.path.exists()



#     cmd_filter_mutect_calls = f' docker run 
#     gatk FilterMutectCalls \
#    -R reference.fasta \
#    -V somatic.vcf.gz \
#    --contamination-table contamination.table \
#    --tumor-segmentation segments.tsv \
#    -O filtered.vcf.gz'





 #As a module

if __name__ == "__main__":
    GATK_module()
    Create_sequence_dictionary()
    Get_pileup_summaries()
#!/usr/bin/env python

import sys
import os
import subprocess
import glob
from os.path import join
import json
import csv
import logging


#Define the function trimming_fastqs

def trimming_fastqs(input_dir, output_dir):

    paired_dictionary = {}

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)


    # Obtain all fastq.gz files and it's path

    fastq_files = glob.glob(f"{input_dir}/*.fastq.gz")
    #sys.exit()

    #Obtain the dictionary of the paired runs for each sample

    sample_list = []
   
    for fastq in fastq_files:
        fastq_name = os.path.basename(fastq)
        sample_name =  fastq_name.split('_')[0]
        if sample_name not in sample_list:
            sample_list.append(sample_name)

        if not sample_name in paired_dictionary:
            paired_dictionary[sample_name] = {}
        if "R1" in fastq:
            paired_dictionary[sample_name]['R1'] = fastq
        if "R2" in fastq:
            paired_dictionary[sample_name]['R2'] = fastq



    for sample_name in paired_dictionary:

        #Create a directory for each sample

        sample_dir = os.path.join(output_dir, sample_name)
        if not os.path.exists(sample_dir):
            sample_dir_trimming = sample_dir +'/trimming'
            os.makedirs(sample_dir_trimming)

        if os.path.exists(sample_dir) and not os.path.exists(f'{sample_dir}/trimming'):
            sample_dir_trimming = sample_dir +'/trimming'
            os.makedirs(sample_dir_trimming)

            

        #fill paired_dictionary 

        if not "R1" in paired_dictionary[sample_name]:
                continue
        if not "R2" in paired_dictionary[sample_name]:
                continue
            
        # Obtain fq1 and fq2 complete path
        #trimmed = 
        
        fq1 = paired_dictionary[sample_name]['R1'] 
        fq2 = paired_dictionary[sample_name]['R2'] 
        trimmed_fq1_name = os.path.basename(fq1.replace(".fastq.gz", ".trimmed.fastq.gz"))
        trimmed_fq2_name = os.path.basename(fq2.replace(".fastq.gz", ".trimmed.fastq.gz"))

        
        # + f'{sample_name}/trimming/'
        # fq1 = join(fq1_inicial,trimmed)
        # fq2 = join(fq2_inicial,trimmed)

        #Obtain the name for the trimmed files

        paired_dictionary[sample_name]['TRIMMED_R1'] = sample_dir + '/trimming/' + trimmed_fq1_name
        paired_dictionary[sample_name]['TRIMMED_R2'] = sample_dir + '/trimming/' + trimmed_fq2_name


        # paired_dictionary[sample_name]['TRIMMED_R1'] = fq1.replace(".fastq.gz", ".trimmed.fastq.gz")
        # paired_dictionary[sample_name]['TRIMMED_R2'] = fq2.replace(".fastq.gz", ".trimmed.fastq.gz")
        # trimmed_fq1_name = os.path.basename(fq1.replace(".fastq.gz", ".trimmed.fastq.gz"))
        # trimmed_fq2_name = os.path.basename(fq2.replace(".fastq.gz", ".trimmed.fastq.gz"))
        
        # print(trimmed_fq1_name, trimmed_fq2_name)
        trimmed_fq1_path = sample_dir +'/trimming/'+trimmed_fq1_name 
        trimmed_fq2_path = sample_dir+'/trimming/'+trimmed_fq2_name
        reports_path = sample_dir+'/trimming/'+sample_name
        
        # print(trimmed_fq1_path, trimmed_fq2_path,reports_path)

        #cmd to execute fastp

        if not os.path.exists(trimmed_fq1_path) and not os.path.exists(trimmed_fq2_path):
            cmd = f'fastp -i {fq1}\
            -I {fq2}\
            -o {trimmed_fq1_path}\
            -O {trimmed_fq2_path}\
            -h {reports_path}.fastp.html\
            -j {reports_path}.fastp.json'
            print(cmd)
            p_fastp = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output = p_fastp.stdout.decode("UTF-8")
            error = p_fastp.stderr.decode("UTF-8")

            if p_fastp.returncode != 0:
                msg = f" ERROR: Could not run fastp for sample {sample_name}"
                logging.error(msg)
                logging.error(error)
                sys.exit(1)
         
               

    return paired_dictionary






#Obtain dictionary with fastp's metrics

input_dir_json = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"

def fastp_metrics(input_dir_json):

    summary_qc_dict = {}
    json_directories = os.listdir(input_dir_json)
    # print(json_directories)

    json_files_list = [] 

    #Fill jason_files_list 

    for json_directory in json_directories:
        json_files = glob.glob(f"{input_dir_json}/{json_directory}/trimming/*.fastp.json")
        json_files_list.extend(json_files)
            
    # print(json_files_list )

    for json_files in json_files_list :

        #Obtain name of the json file and it's path

        json_file =  json_files.split('/')[-1]
        sample_name = json_file.split('.')[0]
        # print(json_file)

        json_path = json_files.split(f'{json_file}')[0]
        # print(json_path)

        os.chdir(json_path)

        with open (json_file) as archivo_json:
            archivo_dict = json.load(archivo_json)
            #print(archivo_dict)
            
            # metrics: q20, q30, total_trimmed_bases, %trimmed_bades, reads_with_adapter, %reads_with_adapter, gc_content

            q20 = archivo_dict['summary']['after_filtering']['q20_rate'] * 100
            q30 = archivo_dict['summary']['after_filtering']['q30_rate'] * 100
            bases_before_filtering = archivo_dict['summary']['before_filtering']['total_bases']
            bases_after_filtering = archivo_dict['summary']['after_filtering']['total_bases']
            total_trimmed_bases = bases_before_filtering - bases_after_filtering 
            trimmed_bases_perc = total_trimmed_bases*100/bases_before_filtering
            total_inicial_reads = archivo_dict['summary']['before_filtering']['total_reads']
            reads_with_adapter = archivo_dict['adapter_cutting']['adapter_trimmed_reads']
            reads_with_adapter_perc = reads_with_adapter*100/ archivo_dict['summary']['before_filtering']['total_reads']
            gc_content = archivo_dict['summary']['after_filtering']['gc_content']

            
            if not sample_name in summary_qc_dict:
                summary_qc_dict[sample_name] = {}

            if not "q20_rate" in summary_qc_dict:
                summary_qc_dict[sample_name]['q20_rate'] = q20
            if not "q30_rate" in summary_qc_dict:
                summary_qc_dict[sample_name]['q30_rate'] = q30
            if not "total_bases_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['total_bases_(M)'] = float(bases_before_filtering)/1000000
            if not "total_trimmed_bases_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['total_trimmed_bases_(M)'] = float(total_trimmed_bases)/1000000
            if not "trimmed_bases_perc" in summary_qc_dict:
                summary_qc_dict[sample_name]['trimmed_bases_%'] = trimmed_bases_perc
            if not "inicial_reads_(M)" in summary_qc_dict:
                summary_qc_dict[sample_name]['inicial_reads_(M)'] = float(total_inicial_reads)/1000000
            if not "reads_with_adapter" in summary_qc_dict:
                summary_qc_dict[sample_name]['reads_with_adapter'] = reads_with_adapter
            if not "reads_with_adapter_%" in summary_qc_dict:
                summary_qc_dict[sample_name]['read_with_adapter_%'] = reads_with_adapter_perc
            if not "gc_content_%" in summary_qc_dict:
                summary_qc_dict[sample_name]['gc_content_%'] = float(gc_content)*100
                        

    #print(summary_qc_dict)

    os.chdir(input_dir_json)

    return summary_qc_dict

    with open('summary_qc_dict.csv', "w", newline="") as f:
        w = csv.DictWriter(f, summary_qc_dict.keys())
        w.writeheader()
        w.writerow(summary_qc_dict)

    f.close()
   




def eliminate_duplicates(sorted_bam_files_list, picard_path):

    removed_duplicate_files_list = []

    for sorted_bam_file in sorted_bam_files_list:

        sorted_bam_file_name = sorted_bam_file.split('/')[-1].split('.bam')[0]
        sorted_bam_file_dir = sorted_bam_file.split(f'{sorted_bam_file_name}')[0]
        sample_name = sorted_bam_file_name.split('_')[0]
        removed_duplicates_dir = f'{sorted_bam_file_dir}{sorted_bam_file_name}.removed_duplicates.bam'
        
        if not os.path.exists(removed_duplicates_dir):
            if os.path.exists(sorted_bam_file):

                cmd_1 = f'java -jar {picard_path}  AddOrReplaceReadGroups -I {sorted_bam_file}\
                -O {sorted_bam_file_dir}/{sorted_bam_file_name}.rg.bam\
                -RGID {sample_name}\
                -LB lib1\
                -PL illumina\
                -PU unit1\
                -SM {sample_name} --CREATE_INDEX'
                cmd_2 = f'java -jar {picard_path} MarkDuplicates\
                -I {sorted_bam_file_dir}/{sorted_bam_file_name}.rg.bam\
                -O {removed_duplicates_dir}\
                -M {sorted_bam_file_dir}/{sorted_bam_file_name}.metrics.txt --REMOVE_DUPLICATES True --CREATE_INDEX'

                print(cmd_1)
                p_rg = subprocess.run(cmd_1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_rg = p_rg.stdout.decode("UTF-8")
                error_rg = p_rg.stderr.decode("UTF-8")

                if p_rg.returncode != 0:
                    msg = f"ERROR: Could not run read group for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_rg)
                    sys.exit(1)
                
               
                    
                print(cmd_2)
                p_md = subprocess.run(cmd_2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                

                output_md = p_md.stdout.decode("UTF-8")
                error_md = p_md.stderr.decode("UTF-8")

                if p_md.returncode != 0:
                    msg = f"ERROR: Could not run mark duplicates for sample {sample_name}"
                    logging.error(msg)
                    logging.error(error_md)
                    sys.exit(1)
           
                else:
                    if removed_duplicates_dir not in removed_duplicate_files_list:
                        removed_duplicate_files_list.append(removed_duplicates_dir)

            if not os.path.exists(sorted_bam_file):
                print(f'Failed to remove duplicates for sample {sample_name}')

        if os.path.exists(removed_duplicates_dir) and not removed_duplicates_dir in removed_duplicate_files_list: 
            removed_duplicate_files_list.append(removed_duplicates_dir)
        
    print(removed_duplicate_files_list)
        

    return removed_duplicate_files_list




#As a module to import:

if __name__ == "__main__":
    trimming_fastqs()
    fastp_metrics()
    eliminate_duplicates()


# To execute in this script:

#if __name__ == "__main__":

    # input_dir = "/home/gencardio/tfm_eli/Mostres_30_ultimes",
    # output_dir = "/home/gencardio/tfm_eli/Mostres_30_ultimes/RESULTATS"

    # fastq_trimming = trimming_fastqs(input_dir, output_dir)
   


   
    # print(paired_dictionary[sample])
    # print(sample, paired_dictionary[sample]['R1'], paired_dictionary[sample]['R1'])

#     print(mostres_name, mostres)

# print(paired_dictionary)
# sys.exit()

# print(mostres_list)
# print(paired_dictionary)


#Create a directory of each sample 
# os.chdir(files) 
# trimmed_directory_list = []

# for mostra in mostres_list :
#     mostra_directory = join(files,mostra)
#     if not os.path.exists(mostra_directory):
#         os.mkdir(mostra_directory)
#     if not mostra_directory in trimmed_directory_list:
#         trimmed_directory_list.append(mostra_directory)
# print(trimmed_directory_list)


# for trimmed_directory in trimmed_directory_list:
#     file_name = trimmed_directory.split('/')[-1]
#     # print(file_name)
#     if file_name in paired_dictionary:
#         names = paired_dictionary.get(file_name)
#         # print(names)
#         file_name_R1 = names.get('R1')
#         file_name_R2 = names.get('R2')
#         html_name = file_name_R1.split('_')[0]
#         # print(file_name_R1)
#         # print(file_name_R2)
    
      


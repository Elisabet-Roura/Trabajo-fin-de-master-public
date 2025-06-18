#!/usr/bin/env python

import sys
import os
import subprocess
import glob
import logging


#Download reference genome and index it

# genome_version = 'hg38'
# download_genome_dir = f"/home/gencardio/tfm_eli/PIPELINE/{genome_version}"


def download_reference_genome(download_genome_dir,genome_version):
    
    #Create the directory to downlad the genome

    if not os.path.isdir(download_genome_dir):
        os.mkdir(download_genome_dir)
    
    os.chdir(download_genome_dir)
    
    if not os.path.exists(f'{download_genome_dir}/{genome_version}.fa.gz'):

        #Download the genome
        download = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{genome_version}/bigZips/{genome_version}.fa.gz'

        #Unzip de genome
        unzip = f'gunzip {genome_version}.fa.gz'

        #Index the genome
        index = f'bwa index {genome_version}.fa'

        subprocess.run(download, shell = True)
        subprocess.run(unzip, shell = True)
        subprocess.run(index, shell = True)
        print(download, unzip, index)
    
    genome_path = f'{download_genome_dir}/{genome_version}.fa'

    
    return genome_path

# genome_version = 'hg38'
# chromosome_name = 'chrM'
# download_genome_dir = f"/home/gencardio/tfm_eli/PIPELINE/{genome_version}/{chromosome_name}"

def download_reference_chromosome(download_genome_dir,genome_version,chromosome_name):
    
    #Create the directory to downlad the genome

    if not os.path.isdir(download_genome_dir):
        os.mkdir(download_genome_dir)
    
    os.chdir(download_genome_dir)
    
    if not os.path.exists(f'{download_genome_dir}/{chromosome_name}.fa.gz'):

        #Download the genome
        download = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{genome_version}/chromosomes/{chromosome_name}.fa.gz'

        #Unzip de genome
        unzip = f'gunzip {chromosome_name}.fa.gz'

        #Index the genome
        index = f'bwa index {chromosome_name}.fa'

        subprocess.run(download, shell = True)
        subprocess.run(unzip, shell = True)
        subprocess.run(index, shell = True)
        print(download, unzip, index)
    
    chrM_path = f'{download_genome_dir}/{chromosome_name}.fa'

    
    return chrM_path

#Funció per obtenir nom dels fixers sorted.bam

def sorted_bam_files(input_fastq, genome_reference):

    final_sorted_bam_list = []

    fastqs_dir = os.listdir(input_fastq)
    print(fastqs_dir)

    for sample_name in fastqs_dir:
        if sample_name.startswith("RB"):

            sample_dir = os.path.join(input_fastq, sample_name)
            sample_dir = sample_dir +'/mapping/'+ genome_reference

            final_file_directory = f'{sample_dir}/{sample_name}_{genome_reference}.sorted.bam'
            final_sorted_bam_list.append(final_file_directory)
    
    print(final_sorted_bam_list)

    return final_sorted_bam_list


#Funció per mapejar cada parella de fastq i generar un fitxer sam

def mapping_fastqs(input_fastq, paired_dictionary, genome_path, genome_reference):


    sam_files_list = []

    print(genome_path)
    

    fastqs_dir = os.listdir(input_fastq)
    print(fastqs_dir)

    
    for sample_name in fastqs_dir:

        if sample_name.startswith("RB"):

            sample_dir = os.path.join(input_fastq, sample_name)
            sample_dir = sample_dir +'/mapping/'+ genome_reference

            if not os.path.isdir(sample_dir):
                os.makedirs(sample_dir)
        
        sam_path = f'{sample_dir}/{sample_name}_{genome_reference}.sam'
        sorted_bam_path = f'{sample_dir}/{sample_name}_{genome_reference}.sorted.bam'
        
        if not os.path.exists(sorted_bam_path):     

            if sample_name in paired_dictionary:
            
                
                trimmed_R1_file = paired_dictionary[sample_name]['TRIMMED_R1'] 
                trimmed_R2_file = paired_dictionary[sample_name]['TRIMMED_R2'] 

                print(trimmed_R1_file,trimmed_R2_file)
                print(genome_path)

                    
                if not os.path.exists(sam_path):

                    cmd = f'bwa mem {genome_path} {trimmed_R1_file} {trimmed_R2_file} -t 4 > {sam_path}'
                    print(cmd)
                    p_map = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    
                    output_map = p_map.stdout.decode("UTF-8")
                    error_map = p_map.stderr.decode("UTF-8")

                    if p_map.returncode != 0:
                        msg = f"ERROR: Could not run mapping command for sample {sample_name}"
                        logging.error(msg)
                        logging.error(error_map)
                        sys.exit(1)

                    else:
                        sam_files_list.append(sam_path)

                else:
                    if os.path.exists(sam_path) and not sam_path in sam_files_list:
                            sam_files_list.append(sam_path) 
        
        if os.path.exists(sorted_bam_path) and not sam_path in sam_files_list:
            sam_files_list.append(sam_path)


    print(sam_files_list)
         
    return sam_files_list



def from_sam_to_bam(sam_files_list, final_sorted_bam_files):

    unsorted_bam_files_list = []
    
    for sorted_bam_file in final_sorted_bam_files:

        sam_path = sorted_bam_file.replace(".sorted.bam", ".sam")

        if sam_path in sam_files_list:

            unsorted_bam_path =  sam_path.replace(".sam",".unsorted.bam")


            if os.path.exists(sorted_bam_file):
                if not unsorted_bam_path in unsorted_bam_files_list:
                    unsorted_bam_files_list.append(unsorted_bam_path)

            if not os.path.exists(sorted_bam_file):
  
                if os.path.exists(sam_path):
                    
                    if not os.path.exists(unsorted_bam_path):
                        if not unsorted_bam_path in unsorted_bam_files_list: 
                            
                            cmd = f'samtools view -bS {sam_path} > {unsorted_bam_path}'
                            print(cmd)
                            p_bam = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            

                            output_bam = p_bam.stdout.decode("UTF-8")
                            error_bam = p_bam.stderr.decode("UTF-8")

                            if p_bam.returncode != 0:
                                msg = f"ERROR: Could not convert SAM to BAM for sample {sample_name}"
                                logging.error(msg)
                                logging.error(error_bam)
                                sys.exit(1)
                            else:
                                unsorted_bam_files_list.append(unsorted_bam_path)

            
                if os.path.exists(unsorted_bam_path):
                    if not unsorted_bam_path in unsorted_bam_files_list:
                        unsorted_bam_files_list.append(unsorted_bam_path)


#   unsorted_bam_files_list = []
    
#     for sorted_bam_file in sorted_bam_files:

#         for sam_file in sam_files_list:

#             sam_file_name = sam_file.split('/')[-1]
#             sam_directory = sam_file.split(f'{sam_file_name}')[0]
#             # print(sam_directory, sam_file_name)

#             unsorted_bam_file_name = sam_file_name.replace(".sam",".unsorted.bam")
#             bam_directory = sam_directory + unsorted_bam_file_name


#             if os.path.exists(sorted_bam_bai_file):
#                 if not bam_directory in unsorted_bam_files_list:
#                     unsorted_bam_files_list.append(f'{bam_directory}')

#             if not os.path.exists(sorted_bam_bai_file):
#                     # print(unsorted_bam_file_name, bam_directory)
#                 if os.path.exists(f'{sam_directory}/{sam_file_name}'):
                    
#                     if not os.path.exists(bam_directory):
#                         if not bam_directory in unsorted_bam_files_list: 
#                             os.chdir(sam_directory)
#                             cmd = f'samtools view -bS {sam_file_name} > {unsorted_bam_file_name}'
#                             subprocess.run(cmd, shell = True)
#                             print(cmd)
#                             unsorted_bam_files_list.append(f'{bam_directory}')

            
#                 if os.path.exists(bam_directory):
#                     if not bam_directory in unsorted_bam_files_list:
#                         unsorted_bam_files_list.append(f'{bam_directory}')


#   unsorted_bam_files_list = []
    
#     for sorted_bam_bai_file in sorted_bam_bai_files:

#         for sam_file in sam_files_list:

#             sam_file_name = sam_file.split('/')[-1]
#             sam_directory = sam_file.split(f'{sam_file_name}')[0]
#             # print(sam_directory, sam_file_name)

#             unsorted_bam_file_name = sam_file_name.replace(".sam",".unsorted.bam")
#             bam_directory = sam_directory + unsorted_bam_file_name


            
#             if os.path.exists(sorted_bam_bai_file) and not {bam_directory} in unsorted_bam_files_list:
#                 unsorted_bam_files_list.append(f'{bam_directory}')

#             if not os.path.exists(sorted_bam_bai_file):
#                     # print(unsorted_bam_file_name, bam_directory)
#                 if os.path.exists(f'{sam_directory}/{sam_file_name}'):
                    
#                     if not os.path.exists(bam_directory):
#                         if not {bam_directory} in unsorted_bam_files_list: 
#                             os.chdir(sam_directory)
#                             cmd = f'samtools view -bS {sam_file_name} > {unsorted_bam_file_name}'
#                             subprocess.run(cmd, shell = True)
#                             print(cmd)
#                             unsorted_bam_files_list.append(f'{bam_directory}')

            
#                 if os.path.exists(bam_directory):
#                     if not {bam_directory} in unsorted_bam_files_list:
#                         unsorted_bam_files_list.append(f'{bam_directory}')
                        
    print(unsorted_bam_files_list)

    return unsorted_bam_files_list
        


 #Funció que ordeni les coordenades del bam (unsorted) i generi un nou bam (sorted) 


def from_unsorted_to_sorted_bam(unsorted_bam_files_list):

    sorted_bam_files_list = []

    for unsorted_bam_file in unsorted_bam_files_list:

        sorted_bam_file = unsorted_bam_file.replace(".unsorted.bam", ".sorted.bam")



        if not os.path.exists(sorted_bam_file): 

            cmd = f'samtools sort -T TEST {unsorted_bam_file} -o {sorted_bam_file}'
            print(cmd)
            p_sort = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_sort = p_sort.stdout.decode("UTF-8")
            error_sort = p_sort.stderr.decode("UTF-8")

            if p_sort.returncode != 0:
                msg = f"ERROR: Could not sort BAM file for sample {sample_name}"
                logging.error(msg)
                logging.error(error_sort)
                sys.exit(1)
            else:
                sorted_bam_files_list.append(sorted_bam_file)

        if os.path.exists(sorted_bam_file): 
            if not sorted_bam_file in sorted_bam_files_list:
                sorted_bam_files_list.append(sorted_bam_file)
    
        
    
    print(sorted_bam_files_list)

    return sorted_bam_files_list



#Funció per indexar fitxers bam.
def index_bam(removed_duplicate_files_list):

    index_bam_files_list = []

    for removed_duplicates_file in removed_duplicate_files_list:

        index_sorted_bam_file = removed_duplicates_file.replace('.sorted.removed_duplicates.bam', '.sorted.removed_duplicates.bam.bai')

        if not os.path.exists(index_sorted_bam_file):

            cmd = f'samtools index {removed_duplicates_file}'
            print(cmd)
            p_index = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_index = p_index.stdout.decode("UTF-8")
            error_index = p_index.stderr.decode("UTF-8")

            if p_index.returncode != 0:
                msg = f"ERROR: Could not index sorted BAM file for sample {sample_name}"
                logging.error(msg)
                logging.error(error_index)
                sys.exit(1)
            else:
                index_bam_files_list.append(index_sorted_bam_file)
                
        if os.path.exists(index_sorted_bam_file):
            if not index_sorted_bam_file in index_bam_files_list:
                index_bam_files_list.append(index_sorted_bam_file)
    
    print(index_bam_files_list)

    return index_bam_files_list   

 #Eliminate .sam files 

        
def remove_files(index_bam_files_list):

    for index_bam_file in index_bam_files_list:
        if os.path.exists(index_bam_file):

            #Eliminate .sam
            sam_file_dir = index_bam_file.replace('.sorted.removed_duplicates.bam.bai','.sam')

            if not os.path.exists(sam_file_dir):
                print(f'{sam_file_dir} does not exist')
            
            else: 
                if os.path.exists(sam_file_dir):
                    os.remove(sam_file_dir)
                    msg = f'{sam_file_dir} has been deleted'
                    print(msg)

            
            
            #Eliminate .unsorted.bam
            unsorted_bam_dir = index_bam_file.replace('.sorted.removed_duplicates.bam.bai','.unsorted.bam')

            if not os.path.exists(unsorted_bam_dir):
                print(f'{unsorted_bam_dir} does not exist')

            else:
                if os.path.exists(unsorted_bam_dir):
                    os.remove(unsorted_bam_dir)
                    msg = f'{unsorted_bam_dir} has been deleted'
                    print(msg)

        
                  


#As a module

if __name__ == "__main__":

    download_reference_genome()
    download_reference_chromosome()
    mapping_fastqs()
    sorted_bam_files()
    from_sam_to_bam()
    from_unsorted_to_sorted_bam()
    index_bam()
    remove_sam_files()

#To try in this script

# if __name__ == "__main__":

#     genome_version = 'hg38'
#     download_genome_dir = f"/home/gencardio/tfm_eli/PIPELINE/{genome_version}"

#     download_reference_genome(download_genome_dir,genome_version)









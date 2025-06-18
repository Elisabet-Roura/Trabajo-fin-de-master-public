#!/usr/bin/env python

import sys
import pandas as pd
import os
import subprocess
import logging




def percent_duplications(index_bam_files_list, summary_qc_dict):

    for index_bam_file in index_bam_files_list:

        duplication_file = index_bam_file.replace('sorted.removed_duplicates.bam.bai','sorted.metrics.txt')
        sample_name = index_bam_file.split('/')[-1].split('_')[0]

        found = False
        with open(duplication_file, 'r') as file: 
            for line in file:
                if "PERCENT_DUPLICATION" in line:
                    found = True
                    continue
                if found:

                    tmp = line.strip().split('\t')[-2].replace(',','.')

                    if sample_name not in summary_qc_dict:
                        print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                        sys.exit(1)

                    summary_qc_dict[sample_name]['duplications_%'] = float(tmp)*100
                    break

    print(summary_qc_dict)
                    

    return summary_qc_dict


def collect_insert_size_metrics(removed_duplicate_files_list, picard_path, metrics_dictionary):

    for removed_duplicate_file in removed_duplicate_files_list:

        metrics_file_name = removed_duplicate_file.replace('sorted.removed_duplicates.bam','insert_size_metrics.txt')
        pdf_histogram_file_name =  removed_duplicate_file.replace('sorted.removed_duplicates.bam','histogram.pdf')
        sample_name = metrics_file_name.split('/')[-1].split('_')[0]



        if not os.path.exists(metrics_file_name):

            cmd = f'java -jar {picard_path} CollectInsertSizeMetrics \
            -I {removed_duplicate_file}\
            -O {metrics_file_name}\
            -H {pdf_histogram_file_name}' 
            print(cmd)  
            p_insert = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_insert = p_insert.stdout.decode("UTF-8")
            error_metrics = p_insert.stderr.decode("UTF-8")

            if p_insert.returncode != 0:
                msg = f"ERROR: Could not run CollectInsertSizeMetrics for sample {sample_name}"
                logging.error(msg)
                logging.error(error_insert)
                sys.exit(1)

        else:
            if os.path.exists(metrics_file_name):
                print(f'CollectInsertSizeMetrics report already exists for {removed_duplicate_file}.')

        found = False
        with open(metrics_file_name, 'r') as file: 
            for line in file:
                if "MEAN_INSERT_SIZE" in line:
                    found = True
                    continue
                if found:

                    mean_insert_size = line.strip().split('\t')[5].replace(',','.')
                    sd = line.strip().split('\t')[6].replace(',','.')

                    metrics_dictionary[sample_name]['mean_insert_size'] = mean_insert_size
                    metrics_dictionary[sample_name]['sd']= sd

                    if sample_name not in metrics_dictionary:
                        print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                        sys.exit(1)

                    break

    print(metrics_dictionary)

    return metrics_dictionary

def collect_aligment_summary_metrics(removed_duplicate_files_list, picard_path, genome_path,metrics_dictionary):


    for removed_duplicate_file in removed_duplicate_files_list:

        metrics_file_name = removed_duplicate_file.replace('sorted.removed_duplicates.bam','aligment_summary_metrics.txt')
        sample_name = metrics_file_name.split('/')[-1].split('_')[0]

        if not os.path.exists(metrics_file_name):

            cmd = f'java -jar {picard_path} CollectAlignmentSummaryMetrics\
            -R {genome_path}\
            -I {removed_duplicate_file}\
            -O {metrics_file_name}'
            print(cmd)
            p_summary = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_summary = p_summary.stdout.decode("UTF-8")
            error_summary = p_summary.stderr.decode("UTF-8")

            if p_summary.returncode != 0:
                msg = f"ERROR: Could not run CollectAlignmentSummaryMetrics for sample {sample_name}"
                logging.error(msg)
                logging.error(error_summary)
                sys.exit(1)
            else:
                if os.path.exists(metrics_file_name):
                    print(f'CollectAlignmentSummaryMetrics report already exists for {removed_duplicate_file}.')

        found = False
        with open(metrics_file_name, 'r') as file: 
            for line in file:
                if "SECOND_OF_PAIR" in line:
                    found = True
                    continue
                if found:

                    total_reads = line.strip().split('\t')[1].replace(',','.')
                    pf_reads = line.strip().split('\t')[2].replace(',','.')
                    mean_read_lenght = line.strip().split('\t')[15].replace(',','.')

                    inicial_reads = metrics_dictionary[sample_name]['inicial_reads_(M)'] 

                    metrics_dictionary[sample_name]['total_mapped_reads_(M)'] = float(total_reads)/1000000
                    metrics_dictionary[sample_name][' %_mapped'] = (float(total_reads)/1000000)*100/float(inicial_reads)
                    metrics_dictionary[sample_name]['pf_reads_(M)'] = float(pf_reads)/1000000
                    metrics_dictionary[sample_name]['mean_read_lenght'] = mean_read_lenght

                    if sample_name not in metrics_dictionary:
                        print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                        sys.exit(1)

                    break

    print(metrics_dictionary)

    return metrics_dictionary

def Create_sequence_dictionary(picard_path,genome_path):

    sequence_dictionary = genome_path.replace('.fa','.dict')

    if not os.path.exists(sequence_dictionary):
        cmd = f'java -jar {picard_path} CreateSequenceDictionary -R {genome_path} -O {sequence_dictionary}'
        subprocess.run(cmd, shell= True)
        print(cmd)
    
    index_reference_file = genome_path.replace('.fa','.fa.fai')


    if not os.path.exists(index_reference_file):
        cmd = f'samtools faidx {genome_path}'
        print(cmd)

        p_faidx = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output_faidx = p_faidx.stdout.decode("UTF-8")
        error_faidx = p_faidx.stderr.decode("UTF-8")

        if p_faidx.returncode != 0:
            msg = f"ERROR: Could not index reference genome with samtools faidx for {genome_path}"
            logging.error(msg)
            logging.error(error_faidx)
            sys.exit(1)

    return sequence_dictionary

def obtain_interval_list(picard_path, sequence_dictionary, bed_file_path):

 # http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList

    interval_list_path = bed_file_path + '.list.interval_list'

    if not os.path.exists(interval_list_path):

        cmd = f'java -jar {picard_path}  BedToIntervalList -I {bed_file_path} -O {interval_list_path} -SD {sequence_dictionary}'
        print(cmd)
        p_bed = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output_bed = p_bed.stdout.decode("UTF-8")
        error_bed = p_bed.stderr.decode("UTF-8")

        if p_bed.returncode != 0:
            msg = f"ERROR: Could not convert BED to interval list for {bed_file_path}"
            logging.error(msg)
            logging.error(error_bed)
            sys.exit(1)

    return interval_list_path
      

def CollectHsMetrics(removed_duplicate_files_list, picard_path, genome_path, interval_list_path, metrics_dictionary):

    genome_index = genome_path + '.fai'
    if not os.path.exists(genome_index):
        cmd = f'samtools faidx {genome_path}'
        print(cmd)
        subprocess.run(cmd, shell=True)
        

    for removed_duplicate_file in removed_duplicate_files_list:


        metrics_file_name = removed_duplicate_file.replace('sorted.removed_duplicates.bam','collect_Hs_metrics.txt')
        sample_name = metrics_file_name.split('/')[-1].split('_')[0]
 
        if not os.path.exists(metrics_file_name):
            cmd = f'java -jar {picard_path} CollectHsMetrics -I {removed_duplicate_file}\
            -O {metrics_file_name}\
            -R {genome_path}\
            -BI {interval_list_path}\
            -TI {interval_list_path}'
            print(cmd)
            p_hsmetrics = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            

            output_hs = p_hsmetrics.stdout.decode("UTF-8")
            error_hs = p_hsmetrics.stderr.decode("UTF-8")

            if p_hsmetrics.returncode != 0:
                msg = f"ERROR: Could not run CollectHsMetrics for sample {sample_name}"
                logging.error(msg)
                logging.error(error_hs)
                sys.exit(1)
            
        else:
            if os.path.exists(metrics_file_name):
               print(f'CollectHsMetrics report already exists for {removed_duplicate_file}.')


        found = False
        with open(metrics_file_name, 'r') as file: 
            for line in file:
                if "MEAN_TARGET_COVERAGE" in line:
                    found = True
                    continue
                if found:

                    mean_target_coverage = line.strip().split('\t')[33].replace(',','.')
                    PCT_TARGET_BASES_1X = line.strip().split('\t')[45].replace(',','.')
                    PCT_TARGET_BASES_2X = line.strip().split('\t')[46].replace(',','.')
                    PCT_TARGET_BASES_10X = line.strip().split('\t')[47].replace(',','.')
                    PCT_TARGET_BASES_20X = line.strip().split('\t')[48].replace(',','.')
                    PCT_TARGET_BASES_30X = line.strip().split('\t')[49].replace(',','.')
                    PCT_TARGET_BASES_40X = line.strip().split('\t')[50].replace(',','.')
                    PCT_TARGET_BASES_50X = line.strip().split('\t')[51].replace(',','.')
                    PCT_TARGET_BASES_100X = line.strip().split('\t')[52].replace(',','.')


                    metrics_dictionary[sample_name]['mean_target_coverage'] = mean_target_coverage
                    metrics_dictionary[sample_name]['Call_rate_1X'] = PCT_TARGET_BASES_1X
                    metrics_dictionary[sample_name]['Call_rate_2X'] = PCT_TARGET_BASES_2X
                    metrics_dictionary[sample_name]['Call_rate_10X'] = PCT_TARGET_BASES_10X 
                    metrics_dictionary[sample_name]['Call_rate_20X'] = PCT_TARGET_BASES_20X 
                    metrics_dictionary[sample_name]['Call_rate_30X'] = PCT_TARGET_BASES_30X
                    metrics_dictionary[sample_name]['Call_rate_40X'] = PCT_TARGET_BASES_40X
                    metrics_dictionary[sample_name]['Call_rate_50X'] = PCT_TARGET_BASES_50X
                    metrics_dictionary[sample_name]['Call_rate_100X'] = PCT_TARGET_BASES_100X




                    if sample_name not in metrics_dictionary:
                        print(f'ERROR: Sample {sample_name} was not found in the metrics dictionary.')
                        sys.exit(1)

                    break

    print(metrics_dictionary)



    return metrics_dictionary


def obtain_metrics_csv(metrics_dictionary, csv_path):

    metrics_file = csv_path  + '/Metrics.csv'


    df = pd.DataFrame.from_dict(metrics_dictionary, orient='index')
    print(df)
    df.to_csv(metrics_file, index=True, decimal='.')

    print(metrics_file)
    
    # df = pd.DataFrame.from_dict(metrics_dictionary, orient='index')
    # df.to_excel('CollectHsMetrics_output.xlsx', index_label='Sample'







#D'on trec cada metrica:
# collect hs: call rate
# insert size: mean insert size i sd https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics
#collect aligmentsummarymetrics: num total reads, pf reads, MEAN_READ_LENGTH  https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard




if __name__ == "__main__":

    percent_duplications()
    Create_sequence_dictionary()
    CollectHsMetrics()
    collect_insert_size_metrics()
    collect_aligment_summary_metrics()
    obtain_metrics_csv



#EN PRINCIPI NO CAL !!!!!

# #MOSDEPTH
#         mosdepth_file_global= index_bam_file.replace('sorted.removed_duplicates.bam.bai','mosdepth.globa.dist.txt')
#         mosdepth_file_summary= index_bam_file.replace('sorted.removed_duplicates.bam.bai','mosdepth.summary.txt')
        
#         if not os.path.exists(mosdepth_file_global) and not os.path.exists(mosdepth_file_summary):
#             cmd = f' mosdepth {mosdepth_file_global} {mosdepth_file_summary} {index_bam_file}'
#             subprocess.run(cmd, shell = True)
#             print(cmd)
        
#         if os.path.exists(mosdepth_file_global) and os.path.exists(mosdepth_file_summary):
#              print(f'El informe de calidad de mosdepth para {index_bam_file} ya existente')



#   Usage: mosdepth [options] <prefix> <BAM-or-CRAM>

# Arguments:

#   <prefix>       outputs: `{prefix}.mosdepth.global.dist.txt`
#                           `{prefix}.mosdepth.summary.txt`
#                           `{prefix}.per-base.bed.gz` (unless -n/--no-per-base is specified)
#                           `{prefix}.regions.bed.gz` (if --by is specified)
#                           `{prefix}.quantized.bed.gz` (if --quantize is specified)
#                           `{prefix}.thresholds.bed.gz` (if --thresholds is specified)

#mosdepth options prefix bam
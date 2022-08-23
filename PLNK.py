#!/usr/bin/env python3
# Kayla Schimke
# Contributors: Roger Volden

import sys
import os
import glob
import re
import argparse
import numpy as np
import mappy as mp
import time as t

VERSION = 'v1.1.0'

def commandline():
    '''Parses arguments.'''
    parser = argparse.ArgumentParser(description='Processes and quantifies reads from an active ONT sequencing run.',
                                     add_help=True, prefix_chars='-',
                                     usage='PLNK.py -i [fast5_dir] [options] [>QC_data]')
    parser.add_argument('-i','--input', type=str, required=True, help='Directory containing fast5 files.')
    parser.add_argument('-o','--output', type=str, default=os.getcwd(), help='Directory to create and store output. Defaults to current directory.')
    parser.add_argument('-S','--samples', type=str, required=True, help="Absolute path to a tsv file containing sample data \
                                                                        Format: Sample_Name\tSplint_Name\tIllumina_5\tIllumina_7")
    parser.add_argument('-s','--splint', type=str, required=True, help='Absolute path to splint fasta')
    parser.add_argument('-a','--adapter', type=str, required=True, help='Absolute path for adapter fasta.')
    parser.add_argument('-c','--config', type=str, required=True, help='Absolute path to config file.')
    parser.add_argument('-C','--consensus', type=str, required=True, help='Path to directory containing C3POa scripts.')
    parser.add_argument('-t','--targets', type=str, required=True, help='Absolute path to bed file containing sites of interest.')
    parser.add_argument('-r','--reference', type=str, required=True, help='Absolute path to reference genome (fasta).')
    parser.add_argument('-n','--threads', type=int, default=1, help='Number of CPU threads to use when processing. Defaults to 1.')
    parser.add_argument('-g','--cuda', type=str, required=True, help='Number or list of numbers specifying GPUs.')
    parser.add_argument('-T','--timer', action='store_true', help='Save runtime metrics to log file.')
    parser.add_argument('-V','--verbose', action='store_true', help='Print tool log text to standard error, otherwise create log files for each tool.')
    parser.add_argument('-v','--version', action='version', version='PLNK '+VERSION, help='Prints the PLNK version.')

    args = parser.parse_args()
    return args

def parseTargets(bedfile, samplefile):
    splint_list = set()
    sample_list = set()
    on_target = {}
    total_reads = {}
    for line in open(samplefile):
        l = line.strip().split('\t')
        if l[0] != 'Name':
            sample = l[0]
            sample_list.add(sample)
            splint_list.add(l[1])
            on_target[sample] = 0
            total_reads[sample] = 0

    coverage_dict={sample:{} for sample in sample_list}
    for line in open(bedfile):
        # maybe add a thing to check for bed file headers
        l = line.strip().split('\t')
        chr = l[0]
        start = int(l[1])
        end = int(l[2])
        for sample in sample_list:
            for base in range(start,end,1):
                chr_base = chr+'_'+str(base)
                if chr_base not in coverage_dict:
                    coverage_dict[sample][chr_base]=0

    return splint_list, sample_list, total_reads, on_target, coverage_dict

def getFileList(query_path):
    '''Takes a path and returns a list of fast5 files excluding the most recent fast5.'''
    file_list = glob.glob(query_path + '/*.fast5')
    if file_list:
        exclude = max(file_list, key=os.path.getctime)
        file_list.remove(exclude)
        file_list.sort(key=lambda x:os.path.getctime(x))
    return file_list

def main():
    script_start = t.time()

    args = commandline()

    print('Vollmers Lab PLNK software {0}\nC3POa path: {1}\nInput path: {2}\nOutput path: {3}\nCPUs: {4}\nGPU device: {5}'.format(VERSION,args.consensus,args.input,args.output,args.threads,args.cuda))

    query_path = args.input
    if not query_path.endswith('/'):
        query_path += '/'
    if not os.path.exists(query_path):
        sys.stderr.write('Input path nonexistent.\n')
        sys.exit(1)

    output_path = args.output
    if not output_path.endswith('/'):
        output_path += '/'
    if not os.path.exists(output_path):
        if args.verbose: sys.stderr.write('Making output directory: {0}\n'.format(output_path))
        os.mkdir(output_path)

    if args.verbose: sys.stderr.write('Parsing target and sample data: {0}, {1}\n'.format(args.targets,args.samples))
    splint_list, sample_list, total_reads, on_target, coverage_dict = parseTargets(args.targets,args.samples)

    if args.verbose: sys.stderr.write('Initializing aligner with reference: {0}\n'.format(args.reference))
    mm_align = mp.Aligner(args.reference, preset='sr')

    if args.verbose: sys.stderr.write('Making tmp directory: {0}\n'.format(query_path + 'tmp/'))
    tmp_path = query_path + 'tmp/'
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)

    if args.verbose and args.timer: sys.stderr.write('Starting live processing while timing\n')
    if args.timer:
        time_points = {}
        time_file = open(output_path+'runtime.log','w')
        time_file.write('Fast5\tCreated\tProcessing Start\tProcessing End\tGuppy\tC3POa\tPostprocessing\tDemultiplexing\tTargeting\n')

    finished = set()
    while True:
        file_list = getFileList(query_path)
        for f5_path in file_list:
            if args.timer: process_start = t.time() - script_start
            f5_filename = f5_path.split('/')[-1]
            f5_basename = f5_filename.split('.')[0]

            if f5_filename in finished:
                continue
            if f5_filename not in os.listdir(tmp_path):
                # symbolic link for the file to not take up extra storage
                if args.verbose: sys.stderr.write('Linking fast5 file: {0}\n'.format(f5_filename))
                os.system('cp {0} {1}'.format(f5_path, tmp_path))

            if args.timer: time_points['guppy'] = t.time() - script_start

            current_path = output_path + f5_basename + '/'
            os.mkdir(current_path)

            if args.verbose:
                sys.stderr.write('Basecalling {0}: {1}\n'.format(f5_filename,current_path))
                os.system('guppy_basecaller -i {0} -s {1} -c dna_r9.4.1_450bps_sup.cfg -x cuda:{2} --compress_fastq 1>&2'.format(tmp_path, current_path, args.cuda))
            else:
                os.system('guppy_basecaller -i {0} -s {1} -c dna_r9.4.1_450bps_sup.cfg -x cuda:{2} --compress_fastq > {3}'.format(tmp_path, current_path, args.cuda, current_path+'guppy_console.log'))
            os.system('rm -f {0}'.format(tmp_path + f5_filename)) # get rid of the link

            if args.timer: time_points['C3POa'] = t.time() - script_start

            consensus_path = current_path + 'Consensus/'
            os.mkdir(consensus_path)

            if args.verbose:
                sys.stderr.write('Consensus Calling {0}: {1}\n'.format(f5_basename,consensus_path))
                os.system('python3 {0} -r {1} -o {2} -c {3} -s {4} -l 500 -d 100 -n {5} -g 400 -co'.format(args.consensus+'/C3POa.py',current_path+'/pass/*.fastq.gz',consensus_path,args.config,args.splint,args.threads))
            else:
                os.system('python3 {0} -r {1} -o {2} -c {3} -s {4} -l 500 -d 100 -n {5} -g 400 -co 2> {6}'.format(args.consensus+'/C3POa.py',current_path+'/pass/*.fastq.gz',consensus_path,args.config,args.splint,args.threads,consensus_path+'consensus_console.log'))

            if args.timer: time_points['Post'] = t.time() - script_start

            for umi in sorted(splint_list):
                if args.verbose:
                    sys.stderr.write('Postprocessing {0} in {1}\n'.format(umi, f5_basename))
                    os.system('python3 {0} -i {1} -o {2} -a {3} -t -n {4} -bt'.format(args.consensus+'/C3POa_postprocessing.py',consensus_path+umi+'/R2C2_Consensus.fasta.gz',consensus_path+umi+'/',args.adapter,args.threads))
                else:
                    os.system('python3 {0} -i {1} -o {2} -a {3} -t -n {4} -bt 2> {5}'.format(args.consensus+'/C3POa_postprocessing.py',consensus_path+umi+'/R2C2_Consensus.fasta.gz',consensus_path+umi+'/',args.adapter,args.threads,consensus_path+'postprocessing_console.log'))

            if args.timer: time_points['Demux'] = t.time() - script_start
            if args.verbose: sys.stderr.write('Demultiplexing {0}\n'.format(f5_basename))
            os.system('python3 {0} {1} {2}'.format('demultiplexBWNsamplesheet.py',args.samples,consensus_path))

            if args.timer: time_points['Targets'] = t.time() - script_start
            for sample in sample_list:
                for name, sequence, quality in mp.fastx_read(consensus_path+'demultiplexed/'+sample+'.fasta'):
                    total_reads[sample] += 1
                    matched = False
                    for hit in mm_align.map(sequence):
                        if hit.is_primary:
                            chr, start, end = hit.ctg, hit.r_st, hit.r_en
                            for base in range(start,end,1):
                                chr_base = chr+'_'+str(base)
                                if chr_base in coverage_dict[sample]:
                                    coverage_dict[sample][chr_base] +=1
                                    matched=True
                    if matched:
                        on_target[sample] += 1

            print('{0}\t{1}'.format(f5_filename,sum(total_reads.values())))
            for sample in sorted(sample_list):
                coverage_list = list(coverage_dict[sample].values())
                percent_total = round((total_reads[sample] / sum(total_reads.values())) * 100, 2)
                percent_on_target = round((on_target[sample] / total_reads[sample]) * 100, 2)
                print('{0}\t{1}\t{2}%\t{3}\t{4}%\t{5}'.format(sample,total_reads[sample],percent_total, \
                                                              on_target[sample],percent_on_target, \
                                                              np.median(coverage_list)))

            if args.timer:
                process_end = t.time() - script_start
                f5_created = os.path.getctime(f5_path) - script_start
                guppy_dur = time_points['C3POa'] - time_points['guppy']
                c3poa_dur = time_points['Post'] - time_points['C3POa']
                post_dur = time_points['Demux'] - time_points['Post']
                demux_dur = time_points['Targets'] - time_points['Demux']
                target_dur = process_end - time_points['Targets']
                time_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(f5_filename,f5_created,process_start,process_end,guppy_dur,c3poa_dur,post_dur,demux_dur,target_dur))

                if args.verbose: sys.stderr.write('Elapsed time processing {0}: {1} seconds\n'.format(f5_filename,process_end-process_start))


            finished.add(f5_filename)
    if args.timer: time_file.close()

if __name__ == '__main__':
    main()

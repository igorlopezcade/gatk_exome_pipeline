import os
from ruffus import *

input_fqs = ['test_data/pe_1.fq', 'test_data/pe_2.fq']

ref_genome = '/resource/genomes/bwa/ucsc.hg19'

def touch_file(fname):
    os.system("touch "+fname)

# for f in input_fqs[0]:
#     touch_file(f)

fq_suffix = '.fq'
fastqc_suffix = '_fastqc'
sam_suffix = '.sam'
bam_suffix = '.bam'
sorted_bam_suffix = '_sorted.bam'
dedupbam_suffix = '_sorted_deduped.bam'
bai_suffix = '.bam.bai'
alignstats_suffix = '.bam.flagstats'
covstats_suffix = '_coverage'

def log_and_run(func):
    def _func(infn, outfn):
        cmd = func(infn, outfn)
        #logger.info("running GATK command: "+cmd)
        print(cmd)
        return os.system(cmd)
    _func.__name__ = func.__name__
    return _func


@transform(input_fqs, suffix(fq_suffix), fastqc_suffix)
@log_and_run
def fastqc(input_fqs, out_dir):
    os.mkdir(out_dir)
    cmd = '/resource/programs/fastqc {} -o {}'.format(input_fqs, out_dir)
    return cmd


@transform(input_fqs, suffix(fq_suffix), sam_suffix)
@log_and_run
def align_with_bwa(input_fqs, sam_file):
    cmd = '/resource/programs/bwa-0.7.12/bin/bwa mem {} {} > {}'.format(ref_genome, input_fqs, sam_file)
    return cmd

@follows(align_with_bwa)
@transform(align_with_bwa, suffix(sam_suffix), bam_suffix)
@log_and_run
def convert_sam_to_bam(sam_file, bam_file):
    cmd = 'samtools view -b  {} > {}'.format(sam_file, bam_file)
    return cmd


@follows(convert_sam_to_bam)
@transform(convert_sam_to_bam, suffix(bam_suffix), sorted_bam_suffix)
def sort_bam(in_bam, out_bam):
    touch_file(out_bam)

@follows(sort_bam)
@transform(sort_bam, suffix(sorted_bam_suffix), dedupbam_suffix)
def mark_duplicate_in_bam(in_bam, out_bam):
    touch_file(out_bam)

@follows(mark_duplicate_in_bam)
@transform(mark_duplicate_in_bam, suffix(dedupbam_suffix), bai_suffix)
def index_bam(in_bam, out_index):
    touch_file(out_index)

@follows(mark_duplicate_in_bam)
@transform(mark_duplicate_in_bam, suffix(bam_suffix), alignstats_suffix )
def get_alignement_stats(in_bam, out_stats):
    touch_file(out_stats)

@follows(mark_duplicate_in_bam)
@transform(mark_duplicate_in_bam, suffix(bam_suffix), covstats_suffix )
def get_coverage_stats(in_bam, out_stats):
    touch_file(out_stats)

pipeline_run([convert_sam_to_bam])
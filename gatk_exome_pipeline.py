
__author__ = 'yuan'

from ConfigParser import ConfigParser
import os
import socket


from ruffus import *
import ruffus.cmdline as cmdline

parser = cmdline.get_argparse(description='Perform exome analysis on alignment files in bam format using GATK.')

parser.add_argument("input_bams", nargs='*')

options = parser.parse_args()

#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)


hostname = socket.gethostname()
sys_cfg = ConfigParser()
sys_cfg.read('{}.sys.cfg'.format(hostname))


input_bams = options.input_bams or ['test_data/sample_a.bam', 'test_data/sample_b.bam']

input_bam_exts = set(map(lambda p:os.path.splitext(p)[1], input_bams))
assert len(input_bam_exts) == 1

bam_ext = list(input_bam_exts)[0]

int_ext = bam_ext + ".interval"
@transform(input_bams, suffix(bam_ext), (bam_ext, int_ext))
def create_realigner_target(bam_fn, out_files):
    int_fn = out_files[1]
    open(int_fn, 'w').write('interval for {}\n'.format(bam_fn))

realignedbam_ext = '.realigned'+bam_ext

@follows(create_realigner_target)
@transform(create_realigner_target, suffix(bam_ext), realignedbam_ext)
def realign_indel(in_files, out_bam):
    in_bam, interval_fn = in_files
    open(out_bam, 'w').write('realigned for {} with {}\n'.format(in_bam, interval_fn))

recal_ext = realignedbam_ext+'.recal_data.grp'
@follows(realign_indel)
@transform(realign_indel, suffix(realignedbam_ext), [realignedbam_ext, recal_ext])
def get_recal_group(in_bam, out_fiels):
    recal_group_fn = out_fiels[1]
    open(recal_group_fn, 'w').write('realigned for {}\n'.format(in_bam))

recalbam_ext = '.preprocessed'+bam_ext
@follows(get_recal_group)
@transform(get_recal_group, suffix(realignedbam_ext), recalbam_ext)
def recalibrate_base(in_files, out_bam):
    in_bam, recal_grp_file = in_files
    open(out_bam, 'w').write('recaled for {} with {}\n'.format(in_bam, recal_grp_file))

vcf_ext = '.vcf'
@follows(recalibrate_base)
@transform(recalibrate_base, suffix(recalbam_ext), vcf_ext)
def call_haplotype(in_bam, out_vcf):
    open(out_vcf, 'w').write('recaled for {}\n'.format(in_bam))

snpvcf_ext = '.snp.vcf'
@follows(call_haplotype)
@transform(call_haplotype, suffix(vcf_ext), snpvcf_ext)
def select_snp_variants(in_vcf, out_vcf):
    open(out_vcf, 'w').write('snp in {}\n'.format(in_vcf))

indelvcf_ext = '.indel.vcf'
@follows(call_haplotype)
@transform(call_haplotype, suffix(vcf_ext), indelvcf_ext)
def select_indel_variants(in_vcf, out_vcf):
    open(out_vcf, 'w').write('indel in {}\n'.format(in_vcf))

hf_snpvcf_ext = '.hardfiltered.snp.vcf'
@follows(select_snp_variants)
@transform(select_snp_variants, suffix(snpvcf_ext), hf_snpvcf_ext)
def hardfilter_snp_variants(in_vcf, out_vcf):
    open(out_vcf, 'w').write('hf snp in {}\n'.format(in_vcf))

hf_indelvcf_ext = '.hardfiltered.indel.vcf'
@follows(select_indel_variants)
@transform(select_indel_variants, suffix(indelvcf_ext), hf_indelvcf_ext)
def hardfilter_indel_variants(in_vcf, out_vcf):
    open(out_vcf, 'w').write('hf indel in {}\n'.format(in_vcf))

merged_vcf_ext = '.hardfiltered.vcf'
@follows(hardfilter_snp_variants, hardfilter_indel_variants)
@collate([hardfilter_snp_variants, hardfilter_indel_variants], regex(r'(.+).hardfiltered.*.vcf'),
         r'\1{}'.format(merged_vcf_ext))
def merge_hf_vcf(in_vcfs, out_vcf):
    open(out_vcf, 'w').write('merged {}\n'.format(','.join(in_vcfs)))


@follows(merge_hf_vcf)
@transform((select_snp_variants, select_indel_variants, hardfilter_indel_variants, hardfilter_snp_variants),
           suffix(vcf_ext), '')
def remove_intermediate_vcfs(in_vcf, out):
    os.remove(in_vcf)

@follows(merge_hf_vcf)
@transform([realign_indel], suffix(bam_ext), '')
def remove_intermediate_bams(in_bam, out):
    os.remove(in_bam)


cmdline.run (options)
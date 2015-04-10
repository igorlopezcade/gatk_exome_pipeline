
__author__ = 'yuan'

from ConfigParser import ConfigParser
import os
import socket


from ruffus import *
import ruffus.cmdline as cmdline

hostname = socket.gethostname()
sys_cfg = ConfigParser()
sys_cfg.read('{}.sys.cfg'.format(hostname))


parser = cmdline.get_argparse(description='Perform exome analysis on alignment files in bam format using GATK.')
parser.add_argument('--gatk', default=sys_cfg.get('program', 'gatk'), help='path to GATK jar file')
parser.add_argument("input_bams", nargs='*')
parser.add_argument('--ref', required=True,
                    help='specify which reference to use. It should be consistant with the reference used in alignment.')
parser.add_argument('--intervals', help="One or more genomic intervals over which to operate. GATK engine parameter.")


options = parser.parse_args()

#  standard python logger which can be synchronised across concurrent Ruffus tasks
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)


ngs_cfg = ConfigParser()
ngs_cfg.read('{}.ngs.{}.cfg'.format(hostname, options.ref))
gatk_bundle = {}
gatk_bundle['a1000g_indels'] = ngs_cfg.get('gatk_bundle', '1000g_indel')
gatk_bundle['a1000g_snps'] = ngs_cfg.get('gatk_bundle', '1000g_snps')
gatk_bundle['omni'] = ngs_cfg.get('gatk_bundle', '1000g_omni')
gatk_bundle['dbsnp'] = ngs_cfg.get('gatk_bundle', 'dbsnp')
gatk_bundle['hapmap'] = ngs_cfg.get('gatk_bundle', 'hapmap')
gatk_bundle['ref_fasta'] = ngs_cfg.get('ref', 'fasta')


java_cmd = "java -jar {}".format(options.gatk)

input_bams = options.input_bams or ['test_data/sample_a.bam', 'test_data/sample_b.bam']
input_bams = options.input_bams or ['test_data/sample.bam']
if options.intervals is None:
    options.intervals = ngs_cfg.get('gatk_bundle', 'exome_targets')

input_bam_exts = set(map(lambda p:os.path.splitext(p)[1], input_bams))
assert len(input_bam_exts) == 1

bam_ext = list(input_bam_exts)[0]

basecmd = "{java_cmd} -R {bundle[ref_fasta]} -L {opt.intervals}".format(opt=options, bundle=gatk_bundle, java_cmd=java_cmd)

def log_and_run(func):
    def _func(*args, **kwargs):
        cmd = func(*args, **kwargs)
        print(cmd)
        return os.system(cmd)
    return _func

@log_and_run
def test_lr(a):
    return 'echo {!s}'.format(a)


int_ext = bam_ext + ".intervals"
@transform(input_bams, suffix(bam_ext), (bam_ext, int_ext))
def create_realigner_target(bam_fn, out_files):
    int_fn = out_files[1]
    cmd = "{basecmd} -T RealignerTargetCreator -o {int_file} -I {inputbam} " \
          "--known {bundle[a1000g_indels]} ".format(basecmd=basecmd, inputbam=bam_fn, int_file=int_fn, opt=options,
          bundle=gatk_bundle, java_cmd=java_cmd)
    logger.info(cmd)
    os.system(cmd)

realignedbam_ext = '.realigned'+bam_ext



@follows(create_realigner_target)
@transform(create_realigner_target, suffix(bam_ext), realignedbam_ext)
def realign_indel(in_files, out_bam):
    in_bam, interval_fn = in_files
    cmd = "{basecmd} -T IndelRealigner -I {inputbam} -o {out_bam} -targetIntervals {interval_fn}".format(basecmd=basecmd,
                                                                                                         inputbam=in_bam,
                                                                                                         interval_fn=interval_fn,
                                                                                                         out_bam=out_bam)
    logger.info(cmd)
    os.system(cmd)

recal_ext = realignedbam_ext+'.recal_data.grp'
@follows(realign_indel)
@transform(realign_indel, suffix(realignedbam_ext), [realignedbam_ext, recal_ext])
@log_and_run
def get_recal_group(in_bam, out_fiels):
    recal_group_fn = out_fiels[1]
    cmd = "{basecmd} -T BaseRecalibrator -I {in_bam} -o {out_fn} -knownSites {bundle[dbsnp]}".format(basecmd=basecmd,
                                                                                                    in_bam=in_bam,
                                                                                                    out_fn=recal_group_fn,
                                                                                                    bundle=gatk_bundle)
    return cmd

recalbam_ext = '.preprocessed'+bam_ext
@follows(get_recal_group)
@transform(get_recal_group, suffix(realignedbam_ext), recalbam_ext)
def recalibrate_base(in_files, out_bam):
    in_bam, recal_grp_file = in_files
    cmd = "{basecmd} -T PrintReads -I {in_bam} -o {out_bam} -BQSR {recal_grp_file}".format(basecmd=basecmd, **locals())
    logger.info(cmd)
    return os.system(cmd)


vcf_ext = '.vcf'
@follows(recalibrate_base)
@transform(recalibrate_base, suffix(recalbam_ext), vcf_ext)
def call_haplotype(in_bam, out_vcf):
    cmd = "{basecmd} -T HaplotypeCaller -I {in_bam} -o {out_vcf} --dbsnp {gatk_bundle[dbsnp]} -stand_emit_conf 10.0 " \
          "-stand_call_conf 30.0".format(basecmd=basecmd, gatk_bundle=gatk_bundle, **locals())
    logger.info(cmd)
    return os.system(cmd)

# snpvcf_ext = '.snp.vcf'
# @follows(call_haplotype)
# @transform(call_haplotype, suffix(vcf_ext), snpvcf_ext)
# def select_snp_variants(in_vcf, out_vcf):
#
#
# indelvcf_ext = '.indel.vcf'
# @follows(call_haplotype)
# @transform(call_haplotype, suffix(vcf_ext), indelvcf_ext)
# def select_indel_variants(in_vcf, out_vcf):
#     open(out_vcf, 'w').write('indel in {}\n'.format(in_vcf))

# hf_snpvcf_ext = '.hardfiltered.snp.vcf'
# @follows(select_snp_variants)
# @transform(select_snp_variants, suffix(snpvcf_ext), hf_snpvcf_ext)
# def hardfilter_snp_variants(in_vcf, out_vcf):
#     open(out_vcf, 'w').write('hf snp in {}\n'.format(in_vcf))
#
# hf_indelvcf_ext = '.hardfiltered.indel.vcf'
# @follows(select_indel_variants)
# @transform(select_indel_variants, suffix(indelvcf_ext), hf_indelvcf_ext)
# def hardfilter_indel_variants(in_vcf, out_vcf):
#     open(out_vcf, 'w').write('hf indel in {}\n'.format(in_vcf))
#
# merged_vcf_ext = '.hardfiltered.vcf'
# @follows(hardfilter_snp_variants, hardfilter_indel_variants)
# @collate([hardfilter_snp_variants, hardfilter_indel_variants], regex(r'(.+).hardfiltered.*.vcf'),
#          r'\1{}'.format(merged_vcf_ext))
# def merge_hf_vcf(in_vcfs, out_vcf):
#     open(out_vcf, 'w').write('merged {}\n'.format(','.join(in_vcfs)))
#
#
# @follows(merge_hf_vcf)
# @transform((select_snp_variants, select_indel_variants, hardfilter_indel_variants, hardfilter_snp_variants),
#            suffix(vcf_ext), '')
# def remove_intermediate_vcfs(in_vcf, out):
#     os.remove(in_vcf)
#
# @follows(merge_hf_vcf)
# @transform([realign_indel], suffix(bam_ext), '')
# def remove_intermediate_bams(in_bam, out):
#     os.remove(in_bam)


cmdline.run(options)
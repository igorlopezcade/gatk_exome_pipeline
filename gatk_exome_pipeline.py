__author__ = 'yuan'

from ConfigParser import ConfigParser
import os
import socket


from ruffus import *
import ruffus.cmdline as cmdline

def zero_file(filename):
    if os.path.exists(filename):
        # save the current time of the file
        time_info = os.stat(filename)
        try:
            f = open(filename, 'w')
        except IOError:
            pass
        else:
            f.truncate(0)
            f.close()
            # change the time of the file back to what it was
            os.utime(filename, (time_info.st_atime, time_info.st_mtime))

hostname = socket.gethostname()
sys_cfg = ConfigParser()
sys_cfg.read('{}.sys.cfg'.format(hostname))


parser = cmdline.get_argparse(description='Perform exome analysis on alignment files in bam format using GATK.')
parser.add_argument('--gatk', default=sys_cfg.get('program', 'gatk'), help='path to GATK jar file')
parser.add_argument("input_bams", nargs='*')
parser.add_argument('--ref', required=True,
                    help='specify which reference to use. It should be consistent with the reference used in alignment.')
parser.add_argument('--intervals', help="One or more genomic intervals over which to operate. GATK engine parameter.")


options = parser.parse_args("--verbose 4 --ref b37 --intervals test_data/sample_target.intervals "
                            "-T remove_realign_interval -T remove_realigned_bam -T remove_read_group_file "
                            "-T remove_intermediate_vcfs".split())

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



int_ext = bam_ext + ".intervals"
@transform(input_bams, suffix(bam_ext), int_ext)
def create_realigner_target(bam_fn, int_fn):
    cmd = basecmd + " -T RealignerTargetCreator -o {int_fn} -I {bam_fn} " \
          "--known {bundle[a1000g_indels]} ".format(bundle=gatk_bundle, **locals())
    logger.info(cmd)
    os.system(cmd)

realignedbam_ext = '.realigned'+bam_ext


@follows(create_realigner_target)
@transform(create_realigner_target, suffix(int_ext), add_inputs([r"\1"+bam_ext]), realignedbam_ext)
def realign_indel(in_files, out_bam):
    interval_fn, in_bam = in_files
    cmd = basecmd + " -T IndelRealigner -I {in_bam[0]} -o {out_bam} -targetIntervals {interval_fn}".format(**locals())
    logger.info(cmd)
    os.system(cmd)

recal_ext = realignedbam_ext+'.recal_data.grp'
@follows(realign_indel)
@transform(realign_indel, suffix(realignedbam_ext), recal_ext)
def get_recal_group(in_bam, recal_group_fn):
    cmd = basecmd + " -T BaseRecalibrator -I {in_bam} -o {recal_group_fn} -knownSites {bundle[dbsnp]}".format(
        bundle=gatk_bundle, **locals())
    logger.info(cmd)
    os.system(cmd)


recalbam_ext = '.preprocessed'+bam_ext
@follows(get_recal_group)
@transform(get_recal_group, suffix(recal_ext), add_inputs([r"\1"+realignedbam_ext]),recalbam_ext)
def recalibrate_base(in_files, out_bam):
    recal_grp_file, in_bam = in_files
    cmd = basecmd + " -T PrintReads -I {in_bam[0]} -o {out_bam} -BQSR {recal_grp_file}".format(**locals())
    logger.info(cmd)
    return os.system(cmd)


vcf_ext = '.vcf'
@follows(recalibrate_base)
@transform(recalibrate_base, suffix(recalbam_ext), vcf_ext)
def call_haplotype(in_bam, out_vcf):
    cmd = basecmd + " -T HaplotypeCaller -I {in_bam} -o {out_vcf} --dbsnp {gatk_bundle[dbsnp]} -stand_emit_conf 10.0 " \
          "-stand_call_conf 30.0".format(gatk_bundle=gatk_bundle, **locals())
    logger.info(cmd)
    return os.system(cmd)


def _get_filtrate_sub_cmd(in_vcf, out_vcf, **name_filter):
    return ' '.join([' -T VariantFiltration -V {} -o {}'.format(in_vcf, out_vcf)] +
                    ['-filter "{}" -filterName {}'.format(v, k) for k,v in name_filter.iteritems() ])


LOWQUAL_NAME_FILTER = {'QD_filter':'QD < 5.0', 'LowCoverage': "DP < 10"}

qfvcf_ext = '.qf.vcf'
@follows(call_haplotype)
@transform(call_haplotype, suffix(vcf_ext), qfvcf_ext)
def filtrate_low_qual(in_vcf, out_vcf):
    cmd = basecmd + _get_filtrate_sub_cmd(in_vcf, out_vcf, **LOWQUAL_NAME_FILTER)
    logger.info(cmd)
    return os.system(cmd)


snpvcf_ext = '.snp.vcf'
@follows(filtrate_low_qual)
@transform(filtrate_low_qual, suffix(qfvcf_ext), snpvcf_ext)
def select_snp_variants(in_vcf, out_vcf):
    cmd = "{basecmd} -T SelectVariants -selectType SNP -selectType MNP --variant {in_vcf} -o {out_vcf}".format(
        basecmd=basecmd, **locals())
    logger.info(cmd)
    return os.system(cmd)

indelvcf_ext = '.indel.vcf'
@follows(call_haplotype)
@transform(call_haplotype, suffix(vcf_ext), indelvcf_ext)
def select_indel_variants(in_vcf, out_vcf):
    cmd = "{basecmd} -T SelectVariants -selectType INDEL --variant {in_vcf} -o {out_vcf}".format(
        basecmd=basecmd, **locals())
    logger.info(cmd)
    return os.system(cmd)

SNP_HARD_FILTER = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"


hf_snpvcf_ext = '.hardfiltered.snp.vcf'
@follows(select_snp_variants)
@transform(select_snp_variants, suffix(snpvcf_ext), hf_snpvcf_ext)
def hardfilter_snp_variants(in_vcf, out_vcf):
    cmd = basecmd + _get_filtrate_sub_cmd(in_vcf, out_vcf, snp_hard_filter=SNP_HARD_FILTER)
    logger.info(cmd)
    return os.system(cmd)

INDEL_HARD_FILTER = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"

hf_indelvcf_ext = '.hardfiltered.indel.vcf'
@follows(select_indel_variants)
@transform(select_indel_variants, suffix(indelvcf_ext), hf_indelvcf_ext)
def hardfilter_indel_variants(in_vcf, out_vcf):
    cmd = basecmd + _get_filtrate_sub_cmd(in_vcf, out_vcf, indel_hard_filter=INDEL_HARD_FILTER)
    logger.info(cmd)
    return os.system(cmd)

merged_vcf_ext = '.hardfiltered.vcf'
@follows(hardfilter_snp_variants, hardfilter_indel_variants)
@collate([hardfilter_snp_variants, hardfilter_indel_variants], regex(r'(.+).hardfiltered.*.vcf'),
         r'\1{}'.format(merged_vcf_ext))
def merge_hf_vcf(in_vcfs, out_vcf):
    cmd = ' '.join([basecmd, '-T CombineVariants', '-o '+out_vcf]+['--variant '+f for f in in_vcfs])
    logger.info(cmd)
    return os.system(cmd)

zeroed_ext = '.is_zeroed'

@follows(merge_hf_vcf)
@transform(create_realigner_target, suffix(int_ext), int_ext+zeroed_ext)
def remove_realign_interval(in_fns, out_fn):
    zero_file(in_fns[0])
    open(out_fn, 'w').close()


@follows(merge_hf_vcf)
@transform((filtrate_low_qual, select_snp_variants, select_indel_variants, hardfilter_indel_variants,
            hardfilter_snp_variants),
           suffix(vcf_ext), vcf_ext+zeroed_ext)
def remove_intermediate_vcfs(in_vcf, out):
    zero_file(in_vcf)
    os.remove(in_vcf+'.idx')
    open(out, 'w').close()

@follows(merge_hf_vcf)
@transform(realign_indel, suffix(realignedbam_ext), realignedbam_ext+zeroed_ext)
def remove_realigned_bam(in_fn, out_fn):
    zero_file(in_fn)
    os.remove(in_fn[:-1]+'i')
    open(out_fn, 'w').close()

@follows(merge_hf_vcf)
@transform(get_recal_group, suffix(recal_ext), recal_ext+zeroed_ext)
def remove_read_group_file(in_fn, out_fn):
    zero_file(in_fn[0])
    open(out_fn, 'w').close()

options.history_file = '.gatk_exome_pipeline.ruffus_history.sqlite'

cmdline.run(options, gnu_make_maximal_rebuild_mode=True, checksum_level=1, touch_file_only=True)



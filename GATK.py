__author__ = 'yuan'

gatk_jar = '/resource/programs/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar'

java_cmd = "java -jar {}".format(gatk_jar)

def RealignerTargetCreator(inputbam, interval, ):
    "{java_cmd} -T RealignerTargetCreator -R {ref_fasta} -o {interval} -I {inputbam} " \
    "--known {a1000g_indels} -L {exome_targets}".format(**locals())

## Whole Exome Sequencing Analysis Procedure

---

1. Preprocess input sequence files
   * If the input sequences contain adaptors, use _cutadapt_ or similar tools to trim them off.  
   * Use _FastQC_ to access the quality of input sequences.
    
2. Align input suquences 
   1. Use _bwa_ to align the sequences to reference genome. 
   2. Use _Picard_ to sort and deduplicate the obtained alignment bam files.
   
3. Variant-call with _GATK_
   The preprocessing of bam files and variant-calling is conducted according to [_GATK_ best practices for DNAseq_](https://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq).
   
4. Annotate variants 
   Evaluate the variants with [_Variant Effect Predictor_](http://www.ensembl.org/info/docs/tools/vep/index.html) or similar tools.
   Annotate the variants with allele frequencies in background population.

5. Run coverage analysis
   Use _DepthOfCoverage_ from _GATK_ to analyze coverage in the targeted regions.
    
6. Select and report variants
   Filter the variants according to criteria by requests (typical filters are functions, genomic context, frequencies.)
   Report the selected variants in VCF, BED or EXCEL format.
      
   

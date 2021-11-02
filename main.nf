#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process adapterTrim {
    conda 'bioconda::seqkit fastp'
    memory '2 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq_F
    path fastq_R

    output:
    path 'input_stats.txt'
    path '*R1.qc.fastq.gz', emit: qc_fastq_F
    path '*R2.qc.fastq.gz', emit: qc_fastq_R

    shell:
    '''
    seqkit stats !{fastq_F} -T -j 64 > input_stats.txt
    F=`basename !{fastq_F} _R1_001.fastq.gz`
    fastp -w 32 -q 30 -i "$F"_R1_001.fastq.gz -I "$F"_R2_001.fastq.gz -o "$F"_R1.qc.fastq.gz -O "$F"_R2.qc.fastq.gz \
		-h fastp_report.html

    '''
}

process variantcall {
    conda 'bioconda::bwa samtools lofreq'
    cpus 12
    publishDir params.outdir, mode: 'copy'

    input:
    path qc_fastq_F
    path qc_fastq_R
    path reference

    output:
    path '*'
    path '*.lofreq.final.bam', emit: bam
    path '*_vars.filt.vcf', emit: vcf

    shell:
    '''
    F=`basename !{qc_fastq_F} _R1.qc.fastq.gz`
    bwa index !{reference}
	bwa mem -t !{task.cpus} !{reference} "$F"_R1.qc.fastq.gz "$F"_R2.qc.fastq.gz > "$F".sam
    samtools fixmate -O bam,level=1 -m --threads !{task.cpus} "$F".sam "$F".fixmate.bam
    samtools sort --threads !{task.cpus} -O bam "$F".fixmate.bam > "$F".bam
    samtools markdup --threads !{task.cpus} -S "$F".bam "$F".dedupe.bam
    lofreq viterbi -f !{reference} "$F".dedupe.bam | samtools sort - --threads !{task.cpus} > "$F".lofreq.realign.bam
	lofreq indelqual --dindel -f !{reference} "$F".lofreq.realign.bam | samtools sort - --threads !{task.cpus} > "$F".lofreq.indel.bam
	lofreq alnqual -b "$F".lofreq.indel.bam !{reference} > "$F".lofreq.final.bam
	samtools index "$F".lofreq.final.bam
	lofreq call-parallel --pp-threads !{task.cpus} --force-overwrite --no-default-filter --call-indels -f !{reference} -o "$F"_vars.vcf "$F".lofreq.final.bam
    lofreq filter -i "$F"_vars.vcf -o "$F"_vars.filt.vcf -v 75
    '''
}

process consensus {
    conda 'bioconda::samtools ivar'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(bam), val(bamID)

    output:
    path '*'

    shell:
    '''
    F=$`basename !{bam} .lofreq.final.bam`

    samtools mpileup -A -d 0 -Q0 !{bam} | 
        ivar consensus -q 20 -t 0 -m 10 -n N -p "$F"
    '''

}
process snpEff {
    container 'quay.io/biocontainers/snpeff:5.0--hdfd78af_1'
    cpus 1
    memory '1 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(vcf), val(vcfID)
    
    output:
    path '*.snpEFF.ann.tsv'

    shell:
    '''
	java -Xmx8g -jar /usr/local/share/snpeff-5.0-1/snpEff.jar NC_045512.2 !{vcf} > !{vcfID}.snpEFF.ann.vcf  

    grep -v "^##" !{vcfID}.snpEFF.ann.vcf | \
	tail -n+2 | \
	cut -f8 | \
	sed 's/|/\t/g' | \
	cut -f1-16 | \
	sed '1i INFO\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\tEXON_INTRON_RANK\tHGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tDISTANCE_TO_FEATURE\tERROR' > !{vcfID}.snpEFF.ann.tmp

	grep -v "^##" !{vcfID}.snpEFF.ann.vcf | \
		cut -f1-7 > !{vcfID}.ann.base.vcf

	paste !{vcfID}.ann.base.vcf !{vcfID}.snpEFF.ann.tmp > !{vcfID}.snpEFF.ann.tsv

	rm !{vcfID}.snpEFF.ann.tmp
	rm !{vcfID}.ann.base.vcf 

    '''
}

process split_vcf {
    conda 'pandas openpyxl'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(ann_tsv), val(ann_tsv_name)

    output:
    path '*'
    
    shell:
    template 'info_splitter.py'
}
    
workflow {
    //fastq_F = channel.fromPath('/Users/anajung/Documents/HandleyLab/data/N743_lineage/fastq/N743_I9314_Courtney_B_1_351_IsolateStock1_UDP0173_GTCCACCGCT_GCTAATAGGA_S11_R1_001.fastq.gz')
    //fastq_R = channel.fromPath('/Users/anajung/Documents/HandleyLab/data/N743_lineage/fastq/N743_I9314_Courtney_B_1_351_IsolateStock1_UDP0173_GTCCACCGCT_GCTAATAGGA_S11_R2_001.fastq.gz')
    fastq_F = channel.fromPath( params.for )
    fastq_R = channel.fromPath( params.rev )
    ref = channel.fromPath('/Users/anajung/Documents/HandleyLab/data/N743_lineage/fastq/NC_045512.2.fasta')
    adapterTrim(fastq_F, fastq_R)
    variantcall(adapterTrim.out.qc_fastq_F, adapterTrim.out.qc_fastq_R, ref)
    //vcfdata = channel.fromPath('/Users/anajung/Documents/HandleyLab/scripts/shotgun_sequencing/modules/out/N743_I9314_Courtney_B_1_351_IsolateStock1_UDP0173_GTCCACCGCT_GCTAATAGGA_S11_vars.filt.vcf').map(vcf -> [vcf, vcf.simpleName])
    //bam= channel.fromPath('N743_I9314_Courtney_B_1_351_IsolateStock1_UDP0173_GTCCACCGCT_GCTAATAGGA_S11.lofreq.final.bam').map(bam -> [bam, bam.simpleName])
    consensus(variantcall.out.bam.map(vcf -> [bam, bam.simpleName]))
    snpEff(variantcall.out.vcf.map(vcf -> [vcf, vcf.simpleName]))
    split_vcf(snpEff.out)
}
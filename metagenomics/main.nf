// ./nextflow run nextflow_work/metagenomics.nf
/* 
 * pipeline input parameters 
 */
params.reads = "/scratch/mk5636/data/AP3MB/10_filtered_L001_R{1,2}_001.fastq.gz"
params.transcriptome = "/scratch/work/cgsb/genomes/Public/Vertebrate_mammalian/Homo_sapiens/Ensembl/GRCh38.p10/Homo_sapiens.GRCh38.dna.toplevel.fa"
params.outdir = "/home/mk5636/nextflow_work/metagenomics/out"
params.kraken_db = "/scratch/mk5636/minikraken_20171019_8GB"

println "reads: $params.reads"
println "transcriptome: $params.transcriptome"
println "outdir: $params.outdir"

ref = file(params.transcriptome)
kraken_db_path = file(params.kraken_db)

Channel 
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .set { read_pairs_ch } 

/*
 * bwa-mem align
 */
process align {
    publishDir "${params.outdir}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, file(reads) from read_pairs_ch
     
    output:
    file "${pair_id}_aligned_reads.sam" into aligned_reads_ch
	
	script:
    """
	module load bwa/intel/0.7.15
	bwa mem -M -t ${task.cpus} ${ref} ${reads[0]} ${reads[1]} > ${pair_id}_aligned_reads.sam
    """
}

/*
 * samtools extract
 */
process extract {
    publishDir "${params.outdir}/extracted", mode:'copy'
	
    input:
    file aligned_reads from aligned_reads_ch
     
    output:
    file "${aligned_reads}_unmapped.sam" into unmapped_ch
	
    script:
    """
    module load samtools/intel/1.6
    samtools view -f 12 -F 256 $aligned_reads > ${aligned_reads}_unmapped.sam
    """
}

/*
 * picard sam2fastq
 */
process sam2fastq {
    publishDir "${params.outdir}/sam2fastq", mode:'copy'
	
    input:
    file unmapped_sam from unmapped_ch
	
    output:
    file "${unmapped_sam}.R{1,2}.fastq" into unmapped_fastq_ch, unmapped_fastq_ch2
    file "${unmapped_sam}.merged.fastq" into unmapped_merged_fastq_ch
	
    script:
    """
    module load picard/2.17.11
    java -jar -Xmx2g \$PICARD_JAR SamToFastq I=$unmapped_sam F=${unmapped_sam}.R1.fastq F2=${unmapped_sam}.R2.fastq
    cat ${unmapped_sam}.R1.fastq ${unmapped_sam}.R2.fastq > ${unmapped_sam}.merged.fastq
    """
}

/*
 * spades
 */
process spades {
    publishDir "${params.outdir}", mode:'copy'
	
    input:
    file(reads) from unmapped_fastq_ch

    output:
    file '*' into spades_out_channel
	
    script:
    """
    module load spades/gnu/3.10.1
    spades.py --meta -t 20 -m 300 -o spades -1 ${reads[0]} -2 ${reads[1]}
    """
}

//unmapped_fastq_ch2.subscribe { println it }

/*
 * kraken
 */
process kraken {
    publishDir "${params.outdir}/kraken", mode:'copy'
	
    input:
    file(reads) from unmapped_fastq_ch2
    file kraken_db_path

    output:
    file "kraken_out.txt" into kraken_out_channel, kraken_out_channel2
    file "kraken_unclassified_out.txt" into kraken_unclassified_out_channel
    file "kraken_classified_out.txt" into kraken_classified_out_channel

    script:
    """
    module load kraken/intel/1.0
    kraken \
    --threads ${task.cpus} \
    --fastq-input \
    --output kraken_out.txt \
    --preload \
    --db $kraken_db_path \
    --check-names \
    --unclassified-out kraken_unclassified_out.txt \
    --classified-out kraken_classified_out.txt \
    --paired ${reads[0]} ${reads[1]}    
    """
}

/*
 * kraken translate
 *
process kraken_translate {
    publishDir "${params.outdir}/kraken-translate", mode:'copy'

    input:
    file(kraken_out) from kraken_out_channel
    file kraken_db_path

    output:
    file "kraken_out.mpa.tsv" into kraken_out_mpa_channel

    script:
    """
    module load kraken/intel/1.0
    kraken-translate \
    --db kraken_db_path \
    --mpa-format \
    $kraken_out > kraken_out.mpa.tsv
    """
}

/*
 * krona
 *
process krona{
    publishDir "${params.outdir}/krona", mode:'copy'

    input:
    file(kraken_out) from kraken_out_channel2

    output:
    file "*" into krona_out_channel

    script:
    """
    module load kronatools/2.7
    ktImportTaxonomy \
    -t 3 \
    -s 4 \
    -o krona \
    $kraken_out
    """
}

/*
 * humann2
 *
process humann2{
    publishDir "${params.outdir}/humann2", mode:'copy'

    input:
    file(unmapped_merged_fastq) from unmapped_merged_fastq_ch

    output:
    file '*' into humann_out_channel
    
    script:
    """
    humann2/0.11.1
    humann2 \
    --threads ${task.cpus} \
    --input $unmapped_merged_fastq \
    --output humann2
    """
}

*/

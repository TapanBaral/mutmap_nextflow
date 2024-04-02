




  if (params.sample_file) {
        Channel
            .fromPath(params.sample_file)
            .splitCsv(header: ['sample_id', 'fastq1', 'fastq2', 'bulk1_R1', 'bulk1_R2', 'bulk2_R1', 'bulk2_R2' ], sep: '\t')
            .map{ row-> tuple(row.sample_id, file(row.fastq1), file(row.fastq2), file(row.bulk1_R1), file(row.bulk1_R2), file(row.bulk2_R1), file(row.bulk2_R2)) }
            .set { input_files }

} 

//Reference genome
def genome_fasta_ch            = Channel.value(params.genome_fasta ? file(params.genome_fasta) : [])

process mutmap {
 container  params.dockerfile
 publishDir "${params.output_dir}", mode: 'copy'
 tag "${sample_id}"
 cpu = params.num_threads
  input:
    tuple  val(sample_id), file(fastq1), file(fastq1), file(bulk1_R1), file(bulk1_R2), file(bulk2_R1), file(bulk2_R2)
    path fasta

  output:
    path("mutmap_results/*") , emit: mutmap_results

  script:
  """
 
  mutmap -r ${fasta} \\
        -c ${fastq1},${fastq1} \\
        -b ${bulk1_R1},${bulk1_R2} \\
        -b ${bulk2_R1},${bulk2_R2} \\
        -n ${task.cpus} \\
        --species 'Rice' \\
        -o mutmap_results
  """
}

workflow mutmap_analysis{
    mutmap(input_files,
            genome_fasta_ch
            )
}
workflow{
    mutmap_analysis()
}

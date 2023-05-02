#!/usr/bin/env nextflow
params.input = '/home/ls/rachelcw/projects/rare_disease/data/bams/'

process BAM2JUNC {

    conda '/home/ls/rachelcw/projects/rare_disease/bam2junc.yml'
    // convert bam to junc file and create a txt file of list of junction sample
    input:
    // path to bam file 
    path p

    output:
    file 'test_juncfiles'

    script:
    """
    for bamfile in `ls $p/*.bam`; do
    echo Converting \$bamfile to \$bamfile.junc
    samtools index \$bamfile 
    regtools junctions extract -s -a 8 -m 50 -M 500000 \$bamfile -o \$bamfile.junc
    echo \$bamfile.junc >> test_juncfiles.txt
    done

    """

}

// process CLUSTERING {
//     //using bash script insted docker?
//     container 'leafcutter'
//     input:
//     file 'junc_list.txt'

//     output:
//     file 'results.txt'

//     script:
//     """
//     # Command to analyze data
//     """
// }

// process SPOT {
//     container 'spot'

//     input:
//     file 'preprocessed.txt' from preprocess

//     output:
//     file 'results.txt'

//     script:
//     """
//     # Command to analyze data
//     """
// }

workflow {
    //input_ch = channel.fromPath(params.input)
    bam_path = BAM2JUNC(params.input)
    //results = analyze()

    // output:
    // file(results)
}

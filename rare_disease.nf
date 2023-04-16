params.input = 'data/*'

process BAM2JUNC {
    
    // convert bam to junc file and create a txt file of list of junction sample
    input:
    file(input) from input_ch

    output:
    file 'junc_list.txt'

    script:
    """
    for file in `ls /home/ls/rachelcw/projects/LEAFCUTTER/ccle/junctions/`
    do
    path='/home/ls/rachelcw/projects/LEAFCUTTER/ccle/junctions/'
    fullpath=$path$file
    echo $fullpath >> ccle_juncfiles_20230204.txt
    done
    """

}

process CLUSTERING {
    //using bash script insted docker?
    container 'leafcutter'
    input:
    file 'junc_list.txt'

    output:
    file 'results.txt'

    script:
    """
    # Command to analyze data
    """
}

process SPOT {
    container 'spot'

    input:
    file 'preprocessed.txt' from preprocess

    output:
    file 'results.txt'

    script:
    """
    # Command to analyze data
    """
}

workflow {
    input_ch = channel.fromPath(params.input)

    results = analyze()

    output:
    file(results)
}

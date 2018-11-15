# author: cmdoret
# 2018-11-15

# This is a snakemake Snakefile, written using snakemake 3.13.3

localrules: all, filter_events

################################################################################
# USER CONFIG
# Change these variables to match datafiles or dependencies path
# or just the program name if you have them in your PATH
REF_DIR = 'data/' # Directory containing reference genome
REF_FILE = 'chr01.fa'
HIC_FASTQ = 'data/hic.end{mate}.fastq.gz' # path to fastq files containing Hi-C reads. e.g. file.end{mate} if files are file.end1 and file.end2 
MATES=["1", "2"]
ENZYME = "DpnII" # Restriction enzyme used for the Hi-C experiment.
FILTER_LIBRARY_EVENTS = False
################################################################################
import sys,os
print(sys.version)

# Get longest common prefix between fastq files
SAMPLE_NAME = os.path.commonprefix(expand(HIC_FASTQ, mate=MATES))
BASE_IDX = REF_DIR + 'bowtie/' + REF_FILE
def get_final_reads(f):
    # Returns the path to the reads used to generate the matrix
    if FILTER_LIBRARY_EVENTS:
        return 'data/map.dat.indices.filtered'
    else:
        return 'data/map.dat.indices'


rule all:
    input:
        get_final_reads


rule logging:
    input:
        aligned = 'data/map.dat',
        used = get_final_reads
    log:
        'data/' + SAMPLE_NAME + '_' + REF_FILE + '.log'
    shell:
        """
        reads_aligned=$(wc -l {input.aligned})
        reads_used=$(wc -l {input.used})
        echo "$reads_aligned aligned reads with MQ > 30" >> {log}
        echo "$reads_used reads used to generate matrix" >> {log}
        """


if FILTER_LIBRARY_EVENTS:

    rule filter_events:
        message: "Filter or remove non-chimeric reads. Choose min number of restriction fragments thresholds for -/+ and +/- events."
        input:
            "data/map.dat.indices"
        output:
            "data/map.dat.indices.filtered"
        shell:
            "python python_codes/library_events_original.py {input}"

rule fragment_attribution:
    params:
        enzyme = ENZYME, 
        ref_dir = REF_DIR
    input:
        'data/map.dat'
    output:
        'data/map.dat.indices'
    shell:
        "python python_codes/fragment_attribution.py {input} {params.enzyme} {params.ref_dir}"

rule filter:
    input:
        expand('data/map{mate}.sam.0.sorted', mate=MATES)
    output: 
        'data/map.dat'
    shell:
        "paste {input} | awk '{{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$4,$3,$7,$9,$8}}' > {output}"

rule sort:
    threads: 8
    input:
        "data/map{mate}.sam.0"
    output:
        temp("data/map{mate}.sam.0.sorted")
    shell:
        "sort -k1 -V --parallel={threads} {input} -T data/tmp > {output}"

rule iterative_mapping:
    message: "{threads} threads allocated for iterative alignment of the following files {input}."
    threads: 8
    params:
        index = REF_DIR + 'bowtie/' + REF_FILE,
        base1 = 'data/map' + MATES[0],
        base2 = 'data/map' + MATES[1]
    input:
        fastq = expand(HIC_FASTQ, mate=MATES),
        done_flag = BASE_IDX + '.done'
    output:
        temp('data/map' + MATES[0] + '.sam.0'),
        temp('data/map' + MATES[1] + '.sam.0')
    shell:
        "python python_codes/iterative_alignment2.py {input.fastq} -i {params.index} -p {threads} -o1 {params.base1} -o2 {params.base2}"

rule index_reference:
    threads: 8
    params:
        base_idx = BASE_IDX
    input:
        ref = REF_DIR + REF_FILE
    output:
        # expand(BASE_IDX + ".{ext}", ext=['1.bt2','2.bt2','3.bt2','4.bt2','rev1.bt2','rev2.bt2'])
        done_flag = touch(BASE_IDX + '.done')
    shell:
        "bowtie2-build {input.ref} {params.base_idx}"


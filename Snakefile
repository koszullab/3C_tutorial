# author: cmdoret
# 2018-11-15

# This is a snakemake Snakefile, written using snakemake 3.13.3

localrules: all, filter_events

################################################################################
# USER CONFIG
# Change these variables to match datafiles
SAMPLE = "yGL1635"
REF_FILE = 'data/chr01.fa'
TMPDIR = 'data/tmp/'
OUTDIR = 'data/out/'
HIC_FASTQ = 'data/hic.end{mate}.fastq.gz' # path to fastq files containing Hi-C reads. e.g. file.end{mate} if files are file.end1 and file.end2
MATES=["1", "2"]
ENZYME = "DpnII" # Restriction enzyme used for the Hi-C experiment.
FILTER_LIBRARY_EVENTS = True
################################################################################
import sys,os
from datetime import datetime
print(sys.version)

now = datetime.now().isoformat('_', 'seconds').replace('-','')
LOG_FILE = os.path.join(OUTDIR, 'log', SAMPLE + '_' + now + '.log')


def get_final_reads(f):
    # Returns the path to the reads used to generate the matrix
    ext = '.dat.indices'
    if FILTER_LIBRARY_EVENTS:
        ext +='.filtered'
    return os.path.join(OUTDIR, SAMPLE + ext)


rule all:
    input:
        get_final_reads, LOG_FILE


rule logging:
    input:
        aligned = os.path.join(OUTDIR, '{}.dat'.format(SAMPLE)),
        used = get_final_reads
    output:
        LOG_FILE
    shell:
        """
        reads_aligned=$(wc -l {input.aligned} | awk '{{print $1}}')
        reads_used=$(wc -l {input.used} | awk '{{print $1}}')
        echo "$reads_aligned aligned reads with MQ > 30" >> {output}
        echo "$reads_used reads used to generate matrix" >> {output}
        """


if FILTER_LIBRARY_EVENTS:

    rule filter_events:
        message: "Filter or remove non-chimeric reads. Choose min number of restriction fragments thresholds for -/+ and +/- events."
        input:
            os.path.join(OUTDIR, "{}.dat.indices".format(SAMPLE))
        output:
            os.path.join(OUTDIR, "{}.dat.indices.filtered".format(SAMPLE))
        shell:
            "python python_codes/filter_events.py {input} {output}"

rule fragment_attribution:
    params:
        enz = ENZYME,
        ref = REF_FILE
    input:
        os.path.join(OUTDIR, '{}.dat'.format(SAMPLE))
    output:
        os.path.join(OUTDIR, '{}.dat.indices'.format(SAMPLE))
    shell:
        "python python_codes/fragment_attribution.py -e {params.enz} -r {params.ref} {input} {output}"

rule filter_MAPQ:
    input:
        expand(os.path.join(TMPDIR, 'map{mate}.sam.0.sorted'), mate=MATES)
    output:
        os.path.join(OUTDIR, '{}.dat'.format(SAMPLE))
    shell:
        "paste {input} | awk '{{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$4,$3,$7,$9,$8}}' > {output}"

rule sort:
    threads: 8
    input:
        os.path.join(TMPDIR, "map{mate}.sam.0")
    output:
        temp(os.path.join(TMPDIR, "map{mate}.sam.0.sorted"))
    shell:
        "sort -k1 -V --parallel={threads} {input} -T data/tmp > {output}"

rule iterative_mapping:
    message: "{threads} threads allocated for iterative alignment of the following files {input}."
    threads: 8
    params:
        tmp = TMPDIR,
        ref = REF_FILE,
    input:
        fastq = expand(HIC_FASTQ, mate=MATES)
    output:
        temp(os.path.join(TMPDIR, 'map' + MATES[0] + '.sam.0')),
        temp(os.path.join(TMPDIR, 'map' + MATES[1] + '.sam.0'))
    shell:
        "python python_codes/iterative_alignment2.py -T {params.tmp} -m -r {params.ref} -p {threads} -o1 {output[0]} -o2 {output[1]} {input.fastq}"

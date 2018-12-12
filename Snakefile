# author: cmdoret
# 2018-11-15

# This is a snakemake Snakefile, written using snakemake 3.13.3

localrules: all, filter_events

################################################################################
# USER CONFIG
# Change these variables to match datafiles
SAMPLE = "test_run"
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
# Keep useful summary statistics of the run into a "log" file.
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
# If desired, filter spurious library events (loops and uncuts). This makes
# the matrix more readable by reducing the diagonal intensity.
    rule filter_events:
        message: "Filter or remove non-chimeric reads. Choose min number of "
                 "restriction fragments thresholds for -/+ and +/- events."
        input:
            os.path.join(OUTDIR, "{}.dat.indices".format(SAMPLE))
        output:
            os.path.join(OUTDIR, "{}.dat.indices.filtered".format(SAMPLE))
        shell:
            "python python_codes/filter_events.py -p -i {input} {output}"

rule fragment_attribution:
# Attributes each read to the corresponding restriction fragment.
    params:
        enz = ENZYME,
        ref = REF_FILE
    input:
        os.path.join(OUTDIR, '{}.dat'.format(SAMPLE))
    output:
        os.path.join(OUTDIR, '{}.dat.indices'.format(SAMPLE))
    shell:
        "python python_codes/fragment_attribution.py -e {params.enz} -r {params.ref} {input} {output}"

rule paste_filter:
# Merge sorted reads into Hi-C pairs and remove pairs where at least one read
# map ambiguously.
    threads: 2
    input:
        expand(os.path.join(TMPDIR, 'map{mate}.sorted.sam'), mate=MATES)
    output:
        os.path.join(OUTDIR, '{}.dat'.format(SAMPLE))
    shell:
        """
        paste <(awk -v OFS="\t" '{{if($2 == 0) {{$2 = "+"}}
                                   else {{$2 = "-"}}
                                   print $3, $4, $5, $2}}' {input[0]}) \
              <(awk -v OFS="\t" '{{if($2 == 0) {{$2 = "+"}}
                                   else {{$2 = "-"}}
                                   print $3, $4, $5, $2}}' {input[1]}) \
          | awk '{{if($3 >= 30 && $7 >= 30) print $1, $2, $4, $5, $6, $8}}'> {output}
        """

rule sort_bed:
# Sort reads by chromosome and position and remove potential supplementary
# alignments in case of chimeric reads, to allow merging them into pairs later.
    threads: 8
    params:
        mate = "{mate}"
    input:
        os.path.join(TMPDIR, "map{mate}.sam")
    output:
        os.path.join(TMPDIR, "map{mate}.sorted.sam")
    shell:
        "samtools sort -n -@ {threads} {input} \
           | samtools view -@ {threads} > {output}"

rule iterative_mapping:
# Truncates reads to 20bp, attempts mapping and extend the reads by 20 bp iteratively.
# This allows to keep reads containing a ligation site.
    message: "{threads} threads allocated for iterative alignment of the following files {input}."
    threads: 8
    params:
        tmp = TMPDIR,
        ref = REF_FILE,
    input:
        fastq = HIC_FASTQ
    output:
        os.path.join(TMPDIR, 'map{mate}.sam')
    shell:
        "python python_codes/iterative_alignment2.py -T {params.tmp} -r {params.ref} -p {threads} -o {output} {input.fastq} -l 40"

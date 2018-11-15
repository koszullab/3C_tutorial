# SORT_AND_FILTER.SH : treats iteratively aligned sam.0 files (by the programm iterative_alignment.py)
#Sort the reads and keep only the ones whose both end have a good mapping quality (standard: MQ>30)
#You must be at the root of data/ and src/ folders. Programs must be in the src/ directory and data in the data/ directory.

f=$1

echo '--------------------
Sort according to the read identification.'
# if sort does not have -V option try -d
sort -V -k1 $f'.end1.sam.0' > $f'.end1.sam.0.sorted'
job1=$!
sort -V -k1 $f'.end2.sam.0' > $f'.end2.sam.0.sorted'
job2=$!
wait $job1 $job2

echo '--------------------
Pairing of both mates in a single file.'
paste $f'.end1.sam.0.sorted' $f'.end2.sam.0.sorted' > $f'.merged'

echo '--------------------
Filter: keep pairs of reads that both have a Mapping Quality above 30 (RName, Pos, Flag x 2)'
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$4,$3,$7,$9,$8}' $f'.merged'  > $f'.dat'

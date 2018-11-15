# RUNALL.SH : process HiC fastq files
#to generate .dat.indices.filtered files ready to be used to build contact maps.
#You must be at the root of src/ and data/ folder. Programs must be in the src/ directory and data in the data/ directory.

f=$1 # the common part of the fastq files name. Their suffix must be ".end1.fastq(.gz)" and ".end2.fastq(.gz)"
ref=$2 #the directory containing the reference genome, either 1 fasta file with all chromosomes or a folder with 1 fasta file for each chromosome (in this case the name of the file must be the same as the name of the chromosome)
index=$3 #bowtie index for alignment
length=$4 #length of the reads
enz=$5 #the restriction enzyme
name_ref=$6 #precise the reference used to align the reads (can be the name of the strain)

echo '----------------------------------------
Beginning of program.
Go at the root of the \"Yeast\" Folder
Usage: bash src/runall_yeasv2.sh FILE_NAME (before ".end1.fastq" and ".end2.fastq") <Ref_genome_directory> <index> <length_of_reads> <restriction_enzyme> <name of reference>'

echo "Which step do you want to start at ?
1 aligning
2 sort and filter
3 attribution of restriction fragment
4 filtering
5 ~ cancel and quit ~"
#start=1
read start #the step where the pipeline begins

if [ $start -le 4 ]; then

	#names
	f1=$f'.end1.fastq'
	f2=$f'.end2.fastq'
	falign=$f'_ref'$name_ref
	freport=$falign'_report.txt'

	#decompress if needed
	recompress='no'
	if [ -f $f'.end1.fastq.gz' ]; then
	    if [ $start -eq 1 ]; then
    	        recompress='yes'
    	        echo '--------------------
    	        Unzipping the files.'
    	        pigz -d $f'.end'*'.fastq.gz'
	    fi
	fi

	#step 1
	if [ $start -le 1 ]; then
    	    echo '--------------------
    	    Align the reads.'
    	    python src/iterative_alignment2.py $f1 $f2 -i $index -p4 -s $length -o1 $falign".end1" -o2 $falign".end2"
	fi

	#step 2
	if [ $start -le 2 ]; then
	    echo '--------------------
	    Sort, filter the reads.'
	    bash src/sort_and_filter.sh $falign
            rm $falign.end*.sam.0* $falign'.merged' 2> /dev/null
	fi

	echo "name of dat file:" $falign".dat"
	aligned=`cat $falign".dat" | wc -l`

	#step 3
	if [ $start -le 3 ]; then
	    echo '--------------------
	    Attributes restriction fragment index to every read.'
	    python src/fragment_attribution.py $falign'.dat' $enz $ref
	    #output: file.dat.indices (RName, Pos, Flag, BinnePosition x 2)
	fi

	#step 4
	if [ $start -le 4 ]; then
	    echo '--------------------
	    Filter to remove as much non-chimeric reads as possible.
	    When the 1st plot appears, choose thresholds for -/+ and +/-.
	    They must be number of restriction sites from which there are as many chimeric reads as noise or more.
	    The total amount of each kind of fragment is shown on a second plot.'
	    python src/library_events_original.py $falign'.dat.indices'
            rm $falign'.dat.indices' 2>/dev/null
	    #output: file.dat.indices.filtered
	fi

	informative=`cat $falign.dat.indices.filtered | wc -l`

	#write counts in a report file
	echo "aligned_reads_MQual>30 $aligned" > $freport
	echo "reads_used_for_matrix $informative" >> $freport

	#recompress files if needed
	if [[ $recompress = 'yes' ]]; then
	    echo '--------------------
	    Rezipping files.'
	    pigz $f1
	    pigz $f2
	fi
fi

echo 'End of program.
----------------------------------------'

# Tutorial for 3C/Hi-C data processing


This document presents several scripts and simple commands for the processing, vizualisation and primary analysis of 3C/Hi-C datasets.
These scripts were used during the INSERM workshop [Capturing chromosome conformation: toward a 3D view of genome regulation](http://ateliersinserm.dakini.fr/accueil-66-16.php), May 9-13, 2016 in Paris at the Pasteur Institute and during the last session of [C3BI Training – Introduction to NGS data analysis] (https://c3bi.pasteur.fr/training-introduction-to-ngs-data-analysis ), 28th of June 2016.

For queries or help getting these running, you can send an email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#building-of-the-contacts-map)
* [Normalization of the data](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#normalization-of-the-data)
* [Computation of genomic distance law](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-genomic-distance-law)
* [Computation of correlation matrices](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-correlation-matrices)
* [Directional Index tool to detect TADs](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#directional-index-tool-to-detect-tads)
* [Decomposition into eigen vectors](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#decomposition-into-eigen-vectors)
* [Use of sparse formalism](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#use-of-sparse-formalism)
* [Miscellaneous](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#miscellaneous)

### Dependencies

These scripts will run on Unix-based systems such as Linux or OS X. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems.
For Windows, you may download Python at https://www.python.org/downloads/windows/. In order to run the Unix commands below, you will need something like Cygwin, Windows Subsystems for Linux, or manually install the command line tools with a package manager like Scoop. Then, a few python modules are necessary for diverses operations on arrays and vizualisation.

#### Python (>=3.6)
* Numpy
* Matplotlib
* pysam
* Scipy
* Biopython
* snakemake

These can be readily installed using the supplied ```requirements.txt``` file:

    pip install -Ur requirements.txt

#### External programs

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](http://www.htslib.org/)
## Raw data extraction and alignment
#### Data
A minimal test dataset is provided in the `data` folder. It consists of a reference fasta file and two fastq files containing the pairs of Hi-C reads.

#### Pipeline overview
The whole pipeline is implemented in a [Snakefile](https://snakemake.readthedocs.io/en/stable/) where file path and variables can easily be changed in the "USER CONFIG" section. Each step of the pipeline is described individually below.

You can run the whole pipeline at once using :

```bash
# change the value of -j to the number of cores on your machine
snakemake -j8
```


#### Alignment
The python script `iterative_alignment2.py` performs iterative mapping by first truncating the reads to a given length (20bp by default) and then iteratively extending and remapping the reads that did not map uniquely. This allows to map reads that would otherwise have a ligation site in the middle. Running `python python_codes/iterative_alignment2.py -h` shows the different arguments it can take.

```
usage: iterative_alignment2.py [-h] -r REFERENCE [-p NB_PROCESSORS]
                               [-o OUT_SAM] [-T TEMPDIR] [-m] [-l MIN_LEN]
                               in_fq

positional arguments:
  in_fq                 The fastq file containing Hi-C reads.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Path to the reference genome, in FASTA format.
  -p NB_PROCESSORS, --nb_processors NB_PROCESSORS
                        number of CPUs used for alignment.
  -o OUT_SAM, --out_sam OUT_SAM
                        Path to the output SAM file for the alignment of
                        in_fq.
  -T TEMPDIR, --tempdir TEMPDIR
                        Directory to write temporary files. Defaults to
                        current directory.
  -m, --minimap2        Use minimap2 instead of bowtie for the alignment.
  -l MIN_LEN, --min_len MIN_LEN
                        Minimum length to which reads should be truncated.
```
For example:

```bash
python python_codes/iterative_alignment2.py \
  -r data/chr01.fa
  -o end1.sam data/hic.end1.fastq.gz
```

The script needs to be run separately on both fastq files and each will return a SAM file containing each read once, including those that are not aligned.

The SAM files must then be sorted by read name, and the header remove to facilitate parsing. This can be achieved using the `samtools` utilities.

```bash
samtools sort -n end1.sam | samtools view > end1.sorted.sam
samtools sort -n end2.sam | samtools view > end2.sorted.sam
```

#### Pairing reads
The sam files can then be parsed and combined into a single file using the GNU coreutils suite of tools to match reads of the same pair. We also filter pairs where both reads have a mapping quality above 30.

```bash
paste <(awk -v OFS="\t" '{if($2 == 0) {$2 = "+"}
                           else {$2 = "-"}
                           print $3, $4, $5, $2}' end1.sorted.sam) \
      <(awk -v OFS="\t" '{if($2 == 0) {$2 = "+"}
                           else {$2 = "-"}
                           print $3, $4, $5, $2}' end2.sorted.sam) \
  | awk '{if($3 >= 30 && $7 >= 30)
            print $1, $2, $4, $5, $6, $8}' > pairs.dat

```
The SAM flags 0 and 16 are converted into + and - strands, respectively and the reads are ordered into a convenient table with one pair per line. For each read, the chromosome, position and strand informations are present.


At this stage, you should have a file containing this information:
```
Sc_chr01 54412 - Sc_chr01 54276 +
Sc_chr01 52956 + Sc_chr01 45745 +
Sc_chr01 69448 - Sc_chr01 69297 +
Sc_chr01 79289 - Sc_chr01 79177 +
Sc_chr01 38824 + Sc_chr01 36185 -
Sc_chr01 131720 + Sc_chr01 131874 -
Sc_chr01 195698 - Sc_chr01 195459 +
Sc_chr01 138110 + Sc_chr01 138287 -

```

You can then assign a restriction fragment to each read using the python code [`fragment_attribution.py`](python_codes/fragment_attribution.py):

```
Assigns restriction fragment ID to Hi-C records in a .dat file.
usage: fragment_attribution.py [-h] -e ENZYME -r REFERENCE
                               input_file [output_file]

positional arguments:
  input_file            Path to the input file containing the coordinates of
                        Hi-C interacting pairs (dat format). Use '-' to read
                        from stdin.
  output_file           Path to the output file (dat.indices format). Defaults
                        to stdout.

optional arguments:
  -h, --help            show this help message and exit
  -e ENZYME, --enzyme ENZYME
                        The restriction enzyme used in the Hi-C protocol. E.g.
                        DpnII.
  -r REFERENCE, --reference REFERENCE
                        Path to the reference genome in FASTA format. Can be
                        multiple files separated by commas.

```

For example
```bash
python python_codes/fragment_attribution.py \
  -e DpnII \
  -r data/chr01.fa \
  pairs.dat \
  pairs.dat.indices
```

This will output a file with two additional column (4th and 8th) containing the indices of the restriction fragments on which the reads are located.

```
Sc_chr01	1665	+	4	Sc_chr01	1838	-	4
Sc_chr01	157259	+	437	Sc_chr01	157541	-	437
Sc_chr01	107877	-	290	Sc_chr01	126038	-	343
Sc_chr01	68529	-	184	Sc_chr01	68276	+	184
Sc_chr01	81348	+	215	Sc_chr01	81527	-	215
Sc_chr01	126315	+	344	Sc_chr01	121347	+	326

```

#### Filtering of the data
A check for the proportion of uncrosslinked events (uncuts, loops...) can be carried out at this stage.  
With these information, one can remove uninformative events more cautiously then just removing all events below 10 kb. This filtering is optional and might be necessary when you want to study the structure of chromatin at very short scales like several kb.
We used the home-made python code [`filter_events.py`](python_codes/filter_events.py):

```
usage: filter_events.py [-h] [-i | -t THRESHOLDS THRESHOLDS] [-p]
                        input_file [output_file]

positional arguments:
  input_file            The file containing the coordinates of Hi-C
                        interacting pairs and, the indices of their
                        restriction fragments (.dat.indices format).
  output_file           Path to the output file (filtered dat.indices).
                        Defaults to stdout.

optional arguments:
  -h, --help            show this help message and exit
  -i, --interactive     Interactively shows plots and asks the user for
                        thresholds. Overrides predefined thresholds. Disabled
                        by default.
  -t THRESHOLDS THRESHOLDS, --thresholds THRESHOLDS THRESHOLDS
                        The minimum number of restriction fragments between
                        reads to consider loops and uncut events. Estimated
                        automatically by default. Must be two integers
                        separated by a space (-t <loop> <uncut>).
  -p, --plot_summary    Output a piechart summarizing library composition.

```

For example:
```bash
python filter_events.py pairs.dat.indices pairs.dat.indices.filtered
```

By default, the program will automatically estimate a threshold to remove loop and uncut events. This threshold represents the minimum number of restriction fragments separating reads in a Hi-C pair. This can be done interactively using the `-i` flag and the proportion of each type of library events can also be displayed using the `-p` flag.

You can have a look to the behaviors of the different pairs of reads taking into account the directions of the reads with this kind of plots:

![alt tag](https://github.com/koszullab/3C_analysis_tools/blob/master/pictures/behavior_events_annotated.png)

At this step, close the window, the program will ask you the two thresholds for the uncuts and loops events that you can determine from the previous plot.
Then, the code counts the different types of events, it gives information about the quality of your library preparation.

![alt tag](https://github.com/koszullab/3C_analysis_tools/blob/master/pictures/PieChart_events.png)

A file where non-crosslinked events have been removed is created and name output_alignment_idpt.dat.indices.filtered.

## Building of the contacts map

From the alignment file, you can build a binned contacts map. We used the home-made python code [`Matrice_Creator.py`](python_codes/Matrice_Creator.py): that converts the alignment file information into a dictionary structure and then into an array of the chosen chromosome(s).
There are two arguments to enter: the path of your output alignment file and the size of the bin (in bp).

Ex:
```bash
python Matrice_Creator.py pairs.indices.filtered  1000000
```

You should have a picture with all the chromosomes :

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/MAT_RAW.png)

In this code, you can modify the list of the chromosome(s) you want to display. For example, you can display only the human chromosome 3 (list_chr = ["chr3"]) you shoud have a picture looking like this:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_RAW.png)

## Normalization of the data
We used the normalization procedure called SCN (for Sequential Component Normalization presented in http://www.biomedcentral.com/1471-2164/13/436).
This procedure assumes that every bin should be detected with the same strength. We divide each line of the previous matrice by its sum then each column by its sum. We reiterate this loop several times. It can be shown that after several iterations, the matrice has converged.
Before the iterations, poor interacting bins must be discarded and considered as non detectable bins (due to experimental biases, unmapple sequences etc).
To do that, a threshold is given as an argument to the SCN function so that every bin with a reads number below this value will be replaced by vectors of zeros. To determine this threshold, you can plot the distribution in the number of reads by doing the histogram in a python terminal with [`histo_r.py`](python_codes/histo_r.py):

```python
from pylab import *
from histo_r import *

MATRICE=loadtxt("chr3_RAW.txt");
histo_r(MATRICE.sum(axis=0),100);
```
The function histo_r(V,N) makes a histogram of the vector V in N bins and plots the result as a bar plot (it is an equivalent of the R function hist() ).  

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/histogram_chr3.png)

To normalise the raw contact map, we used the function scn_func coded in the module [scn_human](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/scn_human.py).
```python
import scn_human

mn=scn_human.scn_func(MATRICE,1500);
imshow(mn**0.2,interpolation="none");
colorbar();
savefig('chr3_NORMALISED.png');
np.savetxt('chr3_NORMALISED.txt',MATRICE);
```
The normalised contacts map should look like this:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_NORMALISED.png)

The matrice elements can be viewed as relative frequencies of contacts between different loci (the scores go now from 0 to 1).
It should be noticed that this procedure does not conserve the symetry property of the matrix. To recover the symetry property, you can use a simple line in python:
```python
mn=(mn+mn.T)/2;
```
where mn.T is the transposed matrice of mn.

## Computation of genomic distance law
This plot is important and must be computed at the very first steps of data processing. It reflects the polymer behaviour of the chromatin and thus allows to check the presence or absence of 3D signal.
It consists in computing the mean number of reads in function of the genomic distance separating the two loci.
The computation thus consists in scanning every diagonal of the matrice and taking the average of this sets of elements.
We coded this computation in the function [distance_law_human](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/distance_law_human.py).

```python
from pylab import *
import distance_law_human

# mn: normalised matrice
d=distance_law_human.dist_law(mn)
plt.plot(d,label="Chr 3",linewidth=3.0);
plt.loglog();
plt.xlabel("Genomic distance in bins of 100kb");
plt.ylabel("Contacts Frequency");
plt.legend();
savefig('chr3_distance_law.png');
```
Example of plot obtained for the chromosome 3:

![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_distance_law.png)



## Computation of correlation matrices
The correlation matrice is often computed in Hi-C data analysis. Even if not always present in final publication, it can render  more visible structures that are present in the data especially domains structures (squares in the contacts maps). It consists in looking for correlation in the shared neighbors between each line and column.
It is simply computed by taking the Pearson Coefficient (or another correlation coefficient) between each line and each column of a matrice coded in the function corrcoef from Numpy.

```python
n1=mn.shape[0];
# Computation of genomic distance law matrice:
MAT_DIST =  np.zeros((n1, n1));
for i in range(0,n1) :
    for j in range(0,n1) :
        MAT_DIST[i,j] =  d[abs(j-i)] ;

#imshow(MAT_DIST**0.2);
MAT2=mn/MAT_DIST;
imshow( MAT2**0.2);

# Computation of the correlation matrice:
mc=np.corrcoef(MAT2);
mc[np.isnan(mc)] = 0.0;

mc2=np.corrcoef(mc);
mc2[np.isnan(mc2)] = 0.0;

imshow( mc2,interpolation='none');
colorbar();
savefig('chr3_CORRELATION2.png');
```

Example of plot for the human chromosome 3 with 100 kb bins and 2 iterations of correlation:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_CORRELATION2.png)


## Directional Index tool to detect TADs
This tool is commonly used in Hi-C data analysis. It looks for change in the directionality between "left vector" and "right vector" at a certain loci in the genome. A change could come from the presence of a border between two different compartments in the genome.
It consists in doing a paired T test between "left vector" and "right vector" on each bin along the genome. The size of the "left vector" and "right vector" is put as a parameter and allows to look for domains structures at a specific scale.
We coded this computation in the function [directional_indice](https://github.com/axelcournac/3C_tutorial/blob/master/python_codes/directional_indice.py).

```python
import directional_indice
import matplotlib.gridspec as gridspec

# Computation and plot of Directional Index:
M = np.corrcoef(mn);
n1 = M.shape[0];
scale=20;
DI = directional_indice.directional(M,scale).T;

# PLOT of contacts Map and Directional Index:
gs = gridspec.GridSpec(2, 1, height_ratios=[5,1] );  #  to have several subplots on the same graph
ax0 = plt.subplot(gs[0]);
ax0.imshow(mn**0.2,interpolation="none")
ax0.set_title("Contacts Map");

ax2 = plt.subplot(gs[1], sharex=ax0);
b1=0;
b2=len(DI.T);
borders = list(range(b1,b2));
borders2 = DI.T
borders2 = array(borders2[:,0]);

ax2.set_xlim([b1, b2]);
ax2.set_ylim([-2.0, 2.0]);
ax2.fill_between( borders, 0,borders2, color='red' );
ax2.fill_between( borders, 0,borders2,borders2<0 ,color='green' );
ax2.set_ylabel("DI scale = 2Mb)");
savefig('chr3_DI_bin100kb.png');
```


Example of plot for the human chromosome 3 with 100 kb bins and DI carried out at the scale of 2Mb:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_DI_bin100kb.png)



## Decomposition into eigen vectors
This geometrical transformation allows to decompose the matrice into eigen values and vectors. It is a way to simplify the data or at least to decrease the dimentionality of the mathematical object.
It can be shown that the first eigen vector corresponds to the 2 compartments partition signal of the genome.
It should be kept in mind that this is a first approximation of the 3D structure of a genome and other or sub-compartments can be detected using higher resolution and/or using other tools of compartment detection.

```python
import matplotlib.gridspec as gridspec
import scipy as sc
from scipy import signal
from scipy.linalg import eig

# Computation of eigen vectors:
(V,D) = sc.linalg.eig(mc);
plot(D[:,0]);

# Plot of contacts Map and Eigen vector:
gs = gridspec.GridSpec(2, 1, height_ratios=[5,1] );
ax0 = plt.subplot(gs[0]);
ax0.imshow(mn**0.2,interpolation="none")
ax0.set_title("Contacts Map");

ax2 = plt.subplot(gs[1], sharex=ax0 );
b1=0;
b2=len(D[:,0]);
borders = list(range(b1,b2));
borders2 = D[:,0]
borders2 = array(borders2[:,0]);

ax2.set_xlim([b1, b2]);
#ax2.set_ylim([-2.0, 2.0]);
ax2.fill_between( borders, 0,borders2, color='darkblue');
ax2.fill_between( borders, 0,borders2,borders2<0 ,color='darkblue' );
ax2.set_xlabel("Position along the genome (in bins of 200kb)");
ax2.set_ylabel("Eigen vector");
savefig('chr3_Eigen_bin100kb.png');
```


Example of plot for the human chromosome 3 with 100 kb bins (zoom at the end of the chromosome):
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/chr3_Eigen_bin100kb.png)

## Use of sparse formalism
An alternative and interesting way to mathematically represent the data is the sparse formalism. It is very relevant for matrices in which most of the elements are zero which is ofter the case for human, mouse contacts maps.
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/sparse.png)

In python, functions implemented for the use of sparse formalism can be found in the module scipy.


## Miscellaneous
From the contacts map, you can do a kind of 4C plot of a position of interest by simply select the corresponding line in the matrix:

```python
plot(mn[ int(37000000/BIN),: ] );
xlabel("Postion along the chromosome in bins of 100kb")
ylabel("Normalised Score");
savefig('4Cplot_chr3.png');
```
Example of 4C plot from the normalised contacts map for the position 37000000 bp on the chromosome 3:
![alt tag](https://github.com/axelcournac/3C_analysis_tools/blob/master/pictures/4Cplot_chr3.png)

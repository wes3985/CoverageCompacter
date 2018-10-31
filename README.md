CoverageCompacter
==============

Next generation sequencing (NGS) is an expensive process but its cost could be significantly reduced if only good quality libraries
were sequenced. CoverageCompacter can be used to generate estimates of the coverage of NGS libraries that have been
sequenced at ultra low depth if those libraries were to undergo further sequencing. This is useful for analysing DNA which
is likely to have regions of allelic dropout such as samples that have undergone whole genome amplification (WGA) from single celled
sequencing experiments or other experiments that use low amounts of starting material such as DNA extracted from laser capture 
microdissection experiments. The output from CoverageCompacter can be used to guide which libraries to shortlist for further sequencing
or producing useful summaries of sequencing depth in the BED format.

Principles and applications
===========================
This software was inspired by work by Zhang et al 2015 (see link to original article below). NGS libraries produced from DNA that has 
undergone WGA are likely to contain amplificaton bias, i.e. areas of the genome which are over/under represented. This bias can be 
estimated from low pass sequencing data based on 2 principles as discussed by Zhang et al. (1) Essentially sequencing a library at ultra low sequencing 
can be thought of as 'sampling' a library, we are sequencing regions which are well represented in the overall fragment pool. (2) WGA tends to produce 
amplicons rangeing between 10-50Kb in size. Infering from these 2 concepts; if a library is sequenced at low depth and 1 or more
reads map to a given loci, then it can be asumed that if more sequencing is done then its likely to cover the 10-50kb region surrounding the loci. 

- Original aritcle: https://www.nature.com/articles/ncomms7822
- Pubmed link: https://www.ncbi.nlm.nih.gov/pubmed/25879913

CoverageCompacter uses the BED format output of Samtools depth (see Samtools docs for detals) and compacts this output into smaller loci. CoverageCompacter 
takes the binsize as an argument and will compact areas of coverage into a single loci surounded by areas of no coverage up to half the size of the binSize 
arguement. The reasoning behind this is that up to half an amplicon telomeric and centromeric of a given read are likely be contained within the library 
but not yet sequenced. Where 2 amplicons are separated by less than half a bin, CoverageCompacter will merge the loci into a single loci. See 'Inputs 
and outputs of the software' examples 1-4 for more detail about how this is implemented. The output files can then be more easily compared with established tools 
such as 'bedtools intersect'.

CoverageCompacter can also be used to compress loci with a minimum depth chosen by the user. This can be suplied as an arguement. The 'NoCov' arguement is set to 0 by 
default but if set to 10 (for example) then loci will be identified as having sufficient coverage if they have greater than 10 reads covering a single base. 
The loci will be compressed into regions with above or below the minimum coverage and separated by the binsize. This is useful for identifying regions that may
require more sequencing, or comparing multiple libraries to identify areas of poor coverage across a range of samples. It can also be used to identify regions 
that have excessive coverage or verifying that targetted sequencing experiments have performed appropriately. It can also be used to identify transcript boundaries 
in RNASeq datasets.


Dependencies
============
- Analysis ready BAM files
- Reference genome e.g. hg38
- Samtools (we used samtools v1.3.1)
- Python 3


Installation
============
Directly from git repo:
```
git clone https://github.com/wes3985/CoverageCompacter.git
```
Or using pip:
```
pip install CoverageCompacter
```

Instructions for usage
================
(1) Run samtools depth over your BAM file(s) to generate a depth file. You can use the following command and launch 
each chromosome to a separate core if you are using a cluster.

- bam_file_list 				= a text file containing a list with the full paths to your bam files.
- path_to_ref   				= the full path to the reference genome.
- chr		   					= the chromosome in the same format as your reference.
- depth_file_full_outpath		= the full file path to write the depth data to.

```
samtools depth -a -q 10 -Q 20 --reference $path_to_ref -r $chr -f $bam_file_list > $depth_file_full_outpath
```

(2) Call CoverageCompacter, this can take some time for large genomes and/or many samples, we recommend launching each chromosome to a different core.

- in_depthfile			= the input depth file produced by samtools depth
- outfile				= the full file path to write the depth data to.
- samples				= a string containing the samples separated by commas, if calling from a python shell: "sample1,sample2,sample3"
- CHROM					= The chromosome covered by the infile e.g "chr1", "1". The infile must span only a single chromosome or loci from sequential chromosomes
						  will be merged. This functionality may be improved in future versions depending on demand.
- binSize				= the presumed size of the amplicons, recommend size is 10000 based on Zhang et al (see link above). Setting the binSize to
						  zero also makes this a useful tool to check for transcript boundaries in RNASeq data.
- NoCov					= the minimum depth at which will be regions will be separated into coverage / no coverage

From a UNIX command line or launch script:
```
python CoverageCompacter $in_depthfile $outfile $samples $CHROM $binSize
```

From a python interface:
```
from CoverageCompacter import CoverageCompacter
CoverageCompacter(depth_file, outfile, samples, CHROM, 10000, 0)
```

Inputs and outputs of the software
==================================

INPUT EXAMPLE: This is the output format of samtools depth which is supplied to CoverageCompacter.
col1 = chromosome, col2 = position in chromosome, col3 = depth at position.
See samtools depth docs for more information.

```
chr1	1	0
chr1	2	0
chr1	3	0
chr1	4	0
chr1	5	0
chr1	6	0
chr1	7	0
chr1	8	0
chr1	9	1
chr1	10	1
chr1	11	1
chr1	12	1
chr1	13	1
chr1	14	1
chr1	15	1
```

OUTPUT EXAMPLE

>chr	start	end	size	firstCoveredBase	lastCoveredBase	meanCoverage	NBasesCovered	DepthSum	coverageFraction
>chr1	1	8	8	None	None	0	0	0	0
>chr1	9	15	7	9	15	1	7	7	1

The output file is a headered file in the bed format. Below is a description of each header:

- chr					= chromosome
- start					= start position of the loci
- end					= end postion of the loci
- size					= total size of the loci
- firstCoveredBase		= the first base in the loci with coverage >0
- lastCoveredBase		= the last base in the loci with coverage >0
- meanCoverage			= mean coverage accross the loci
- NBasesCovered	    	= toal number of bases covered in the loci
- DepthSum				= the total sum of depth in the loci, this is useful for identifying regions of high relative coverage
- coverageFraction		= fraction of the loci that has coverage >0

Below is an example input with a description of what happens when different numbers are supplied to the binsize arguement.

INPUT EXAMPLE 1

>chr1	1	0
>chr1	2	0
>chr1	3	1
>chr1	4	1
>chr1	5	1
>chr1	6	0
>chr1	7	0
>chr1	8	0
>chr1	9	1
>chr1	10	1
>chr1	11	0
>chr1	12	0
>chr1	13	1
>chr1	14	1
>chr1	15	0
>chr1	16	1
>chr1	17	1
>chr1	18	1
>chr1	19	0

OUTPUT of CoverageCompacter with binSize set to '0' 
Notes: 
- areas of no coverage are joined into a single loci in the bed format
- areas of continuous coverage are also joined into single loci

>chr	start	end		size	firstCoveredBase	lastCoveredBase	meanCoverage	NBasesCovered	DepthSum	coverageFraction
>chr1	1		2		2		None				None			0.0				0				0			0.0
>chr1	3		5		3		3					5				1.0				3				3			1.0
>chr1	6		8		3		None				None			0.0				0				0			0.0
>chr1	9		10		2		9					10				1.0				2				2			1.0
>chr1	11		12		2		None				None			0.0				0				0			0.0
>chr1	13		14		2		13					14				1.0				2				2			1.0
>chr1	15		15		1		None				None			0.0				0				0			0.0
>chr1	16		18		3		16					18				1.0				3				3			1.0
>chr1	19		19		1		None				None			0.0				0				0			0.0

OUTPUT of CoverageCompacter with binSize set to '1'
Notes:
- The boundary of the bin is moved to 'binSize/2' and rounded up to the nearest whole number - '1' in this case.
- Loci with coverage will include 1 bp with zero coverage on either side. 
- Where a single base has no coverage between 2 loci these adjacent loci will be merged.

>chr	start	end		size	firstCoveredBase	lastCoveredBase	meanCoverage	NBasesCovered	DepthSum	coverageFraction
>chr1	1		1		1		None				None			0.0				0				0			0.0
>chr1	2		6		5		3					5				0.6				3				3			0.6
>chr1	7		7		1		None				None			0.0				0				0			0.0
>chr1	8		11		4		9					10				0.5				2				2			0.5
>chr1	12		19		8		13					18				0.625			5				5			0.625

OUTPUT of CoverageCompacter with binSize set to '2'
Notes:
- The output is the same as supplying binsize of 1 (2/2=1) because a boundary is formed at half the bin size and a single base cannot be cut in half. 
- This example is to illustrate the behavior, in reality we would expect CoverageCompacter to be used either with binSize = 0 or binSize set to much 
larger values such as 10,000

>chr	start	end		size	firstCoveredBase	lastCoveredBase	meanCoverage	NBasesCovered	DepthSum	coverageFraction
>chr1	1		1		1		None				None			0.0				0				0			0.0
>chr1	2		6		5		3					5				0.6				3				3			0.6
>chr1	7		7		1		None				None			0.0				0				0			0.0
>chr1	8		11		4		9					10				0.5				2				2			0.5
>chr1	12		19		8		13					18				0.625			5				5			0.625

OUTPUT of CoverageCompacter with binSize set to '3'
Notes:
- The area between postion 6 and postion 8 has 3 bp with zero coverage, in cases where splits occur that could reasonably fit into 
 both adjacent loci, the loci is split and the extra non-covered base is allocated to the previous loci

>chr	start	end		size	firstCoveredBase	lastCoveredBase	meanCoverage		NBasesCovered	DepthSum	coverageFraction
>chr1	1		7		7		3					5				0.42857142857142855	3				3			0.42857142857142855
>chr1	8		19		12		9					18				0.5833333333333334	7				7			0.5833333333333334


INPUT EXAMPLE 2:

>chr1	1	0
>chr1	2	0
>chr1	3	1
>chr1	4	2
>chr1	5	1
>chr1	6	0
>chr1	7	0
>chr1	8	0
>chr1	9	1
>chr1	10	1
>chr1	11	2
>chr1	12	2
>chr1	13	1
>chr1	14	1
>chr1	15	0
>chr1	16	1
>chr1	17	1
>chr1	18	1
>chr1	19	0

OUTPUT of CoverageCompacter with binSize set to '0' and noCov set to '1'

>chr	start	end		size	firstCoveredBase	lastCoveredBase	meanCoverage	NBasesCovered	DepthSum	coverageFraction
>chr1	1		3		3		None				None			0.0				0				0			0.0
>chr1	4		4		1		4					4				2.0				1				2			1.0
>chr1	5		10		6		None				None			0.0				0				0			0.0
>chr1	11		12		2		11					12				2.0				2				4			1.0
>chr1	13		19		7		None				None			0.0				0				0			0.0

Future Updates
==============
- Add an argument which can be specified by the user to exclude regions of no coverage from the output file.


Support
=======
Please inform us of any issues at 'w.woollard@ucl.ac.uk'

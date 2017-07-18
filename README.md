# Computational pipeline to detect gene loss events

In order to achieve high specificity, this pipeline implements a number of steps that overcome assembly and alignment issues, and address evolutionary exon-intron structure changes in genes that are conserved.

The starting point of the pipeline is a multiple genome alignment. The output is a list of lost genes.

### Quality masking of genome alignments
mafQualityMasking.pl takes a maf file as an input and outputs another maf where every base whose sequence quality is below the threshold is masked to N.

### Cleaning e-lines (apparently deletions or unaligning sequence) in genome alignments that overlap assembly gaps
mafFilterElines.pl takes a maf file as an input and cleans e-lines that overlaps assembly gaps

### Cleaning those chains in the maf that do not represent true orthologs but instead come from paralogs/pseudogenes: 
#### Step 1: filterMafsForWhiteListedChains.pl
This is pre-processing tool that prepares the bedfiles that are needed to run the parent tool. The "PPS" at the end signifies that this  is a pre-preprocessing script

#### Step 2: filterMafsForWhiteListedChains.pl
filterMafsForWhiteListedChains.pl
filterMafsForWhiteListedChains.pl takes a maf file as an input and and removes sequence lines that most likely come from a processed pseudogene or a paralog

### Realigment with CESAR: 
This process involves several steps: starting from determining intron lengths (getIntronTables.pl), finding introns that are deleted (groupExons.pl), combining exons that are separated by deleted intron(s) and re-aligning them with CESAR (RealignmentPipeline_ExonGroups.pl). 

### intronTablesWrapper.pl 
This is a wrapper that runs all the steps mentioned above
Subsequently, jobs file are prepared by intronTableswrapper.pl that run CESAR on the entire set of genes (RealignmentPipeline_Exons.pl). This jobs will need to be pushed to your compute cluster


### Print pairwise alignments in a tabular format. 
The tabular format is a data structure that is used by geneLossPipeline to output mutations in the pairwise alignments (every reference-query species pair). createTabFiles.pl does this

### Detecting mutations in pairwise alignments
This step involves running geneLossPipeline.pl on the tabular files produced in the above step and finds mutations from the pairwise alignments

### Detect gene losses: 
Once the run of geneLossPipeline.pl is complete, geneLossCaller.pl scans through the mutation profile for every species-gene pair and decides whether to call a gene as lost or not.

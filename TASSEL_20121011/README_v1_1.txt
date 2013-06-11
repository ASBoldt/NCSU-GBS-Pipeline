================================================================================
NCSU Fork of TASSEL for Genotyping by Sequencing

Software version: v1.1

Date: 2013-06-10

Legal information:
This software is provided "AS IS".  Use at your own discretion. Use at your own 
discretion. It s released under the GNU Lesser GPL v 2

Contact information:
Ross Whetten (ross_whetten@ncsu.edu)

================================================================================
The documents "Instructions" and "Amazon EC2 script setup" found in the same
directory as this README provide more details and instructions to setting up
and running this program.
=================================================================================

This program / pipeline was built as a tool for streamlining the processing of 
genotyping by sequencing data. It is a fork of the TASSEL software developed and 
made available by the USDA-ARS Lab with Cornell's Institute for Genomic Diversity 
(http://www.maizegenetics.net/).

This software takes two FASTQ files representing forward and reverse sequencing 
of the same flowcell lane and:
1) concatenates the forward and reserve sequencing utilizing the forward 
	sequence identifier
2) logs tag sequences and barcode combinations 
3) generates a tag by taxa output file
4) generates a summary file of the number of forward, reverse, and paired 
	tags at set intervals so that the user can see trends resulting from
	processing the data file

================================================================================
HARDWARE AND SOFTWARE REQUIREMENTS

Hardware: Required (optional)  

OS: 64-bit Linux (Ubuntu)
RAM: 16GB (32+GB)
PROCESSOR: multi-core
STORAGE: A minimum size of space at least 1/4 the combined size of the 
	files you are processing.  Example: Two 40GB compressed files [80GB total]
	, have 20GB free.

Software:
* Oracle Java 1.7 (http://www.java.com)
* git 1.7+ (http://git-scm.com/)
* cd-hit-est for clustering of tags (available in Debian-Med repository and using 
the Amazon EC3 setup script) 

================================================================================
FILE PREPARATIONS

KeyFile Example:
The key file will need to be in the following tab-delimited format:
	Flowcell	Lane	Barcode	Sample	PlateName	Row	Column	 Blank
	C116KACXX	7	ACTCCACG	gbs	GBS1	A	1	Blank
	…

FASTQ Files:
FASTQ files must be named using the following convention
	read#_flowcell_s_lane#_fastq.gz
	Example:  read1_C116JACXX_s_7_fastq.gz
================================================================================
	
INSTALLATION INSTRUCTIONS

To Install:
1. Verify java and git are installed by typing the following from the command prompt
	java -version
	git --version
2. Make a directory where you want to store the program (example: PROGRAM_DIR)
3. Navigate to that directory from the command line
4. Execute the following command within the direcotry:
	git clone git://github.com/conceptstailored/NCSU-GBS-Pipeline.git
5. This should have copied a directory called "NCSU-GBS-Pipeline" in your PROGRAM_DIR
6. Create another directory in your
	 PROGRAM_DIR/NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics called "outputs"
7.  Within PROGRAM_DIR/.../outputs create two directories: 
	/fastq (and copy your fastq files here)
	/FastqPariedEndTagsAndTaxa (where program outputs will be stored)
8. Place your key files in  PROGRAM_DIR/.../outputs

** This is the bare bones approach to running only the paired end plug-in. 
The TASSEL program can do many other things but requires a bit more setup
that can be read about at 
http://www.maizegenetics.net/tassel/docs/TasselPipelineGBS.pdf

--------------------------------------------------------------------------------
DOCUMENTATION:

* Fastq files must be named using the following template
	read#_flowcell_s_lane#_fastq.gz
	
	Example:
	read1_C116JACXX_s_7_fastq.gz

* To run the paired end plugin, execute the following command (all on one line)
	via command line from the PROGRAM_DIR/.../outputs directory.  Optional arguements
	are in () and explained in GETTING STARTED TIPS.

	PROGRAM_DIR/.../maizegenetics/run_gbs_pairedEnd_pipeline.pl -[max_memory] 
	-fork1 -[plugin] -i [input_directory] -k [keyfile1:keyfile2] 
	-e [enzyme1:enzyme2] (-s [maximum_count]) (-c [minimum_detected]) 
	-o [output_directory] -endPlugin -runfork1

	Example:
	PROGRAM_DIR/.../maizegenetics/run_gbs_pairedEnd_pipeline.pl -Xmx30g -fork1 
	-FastqPairedEndToTagCountPlugin -i fastq -k GBS1.key:GBS2.key 
	-e PstI-MspI:MspI-PstI -s 10000 -c 10 -o FastqPairedEndTagsAndTaxa 
	-endPlugin -runfork1
	
* Alternatively, a bash script with the above command, saved and executed in 
	the "outputs" directory reduces typing errors.
	
* The output files you WANT TO KEEP:
	* Barcode_ID_Info.txt = number of barcode combinations detected across all lanes
	* Summary_Stats.txt = printout of the number of forward, reverse, and pairs
		per lane over X interval
	* Tags_Totaled.txt = number of unique tags across all lanes
	* Tags_by_Taxa.txt = number of tags and taxa combination by lane over all lanes

* The plug-in does create large intermediate files that can be deleted once the 
	program completes.  Any files in PROGRAM_DIR/outputs/FastqPariedEndTagsAndTaxa
	that start with "read...." and end with ".txt" can be removed to recover disk
	space.
	
* The plug-in alters the second key file and stores it as [keyfile2]_mod.key.  
	Once processing of all data files is complete, this file can and should be 
	deleted.  If it is left within the directory it will be appended to each time
	the plug-in is run.

--------------------------------------------------------------------------------
GETTING STARTED TIPS:

* Read through the entire README files before attempting to start

* Close any other programs to free system resources

* The program will take several hours to run, and will take longer with each lane
	of data that you process.
** Mandatory flags to set in the calling command

	-km (yes/no) Specify if you want to modify your read2 barcodes.  This flag calls an algorithm that generates permutations of your read2 barcodes by substitution every other base at each position, and deleting one base at each position.  The original barcode is included intact. Use this if your sequence is of poorer quality and you feel that you are missing hits because of a sequencing miscall or error.  WARNING: If your barcodes are used as identifiers and are not of sufficient uniqueness, using this option will create duplicate barcodes and will potential create problems in clearly identifying your sample.

	-wiid (1-5) = where-is-identification = The following are the options that can be used here:
		1 = Barcode in read1 only
		2 = Barcodes in both read2 and read2
		3 = Barcode in read1 and Illumina indexing (CASAVA 1.4-1.7 format) in read2 header
		4 = Barcode in read1 and Illumina indexing (CASAVA 1.8 format) in read2 header
		5 = Barcode in read1 only, sequence generally poor, just grab first 64 bases regardless (this is for diagnostic purposes if your sequence is poor)

	-r1n (yes/no) = Accept missing or uncertain bases in tag area for read 1.
  Determines level of stringency when scanning for barcodes.  Dependent on your barcode keyfile since even if you set this to yes (less stringent) the algorithm still looks for exact matches to your keyfiles, so you would have to add permutations of your barcode with Ns present.

	-r2n  (yes/no)= Accept missing or uncertain bases in tag area for read 2.  Determines level of stringency when scanning for barcodes.  Dependent on your barcode keyfile since even if you set this to yes (less stringent) the algorithm still looks for exact matches to your keyfiles, so you would have to add permutations of your barcode with Ns present.  This can be used in conjunction with the -km flag to allow for 1 N value to be present in your keyfile automatically.  If you want to allow for more than 1 N you will have to make the additions and adjustments to your keyfile.

*  Optional flags that can be set in the calling command
	-s = maximum_count = the default setting is to try and process the 
		entire data set per lane.  This can be extremely memory intensive.  
		If the program crashes due to memory issues, you can look to the last 
		screen output and set this flag to be a maximum number of tags to accept per lane.  
		Using the screen output is only a guide, you will have to play around 
		with this value based on your system's configuration
	
	-c = minimum_detected = the default setting is to accept every tag that is 
		detected.  This may lead to a lot of singlets.  Setting this argument
		can reduce the overall size of your outputs, and can be used to narrow
		your outputs to only a count you are interested in
		
* The keyfile argument flag (-k) needs to have 2 arguments even if you are 
	using the same keyfile for both reads.  The name doesn't matter, so you 
	could enter file1.key:file1.key and that would be accepted.
	
* Enzymes and keyfiles need to be separated by a ":"


================================================================================
Important known bugs, problems, or limitations:

2013-02-21: 
* Software has only been tested with 2 lanes of data (4 read files).
	In theory it should be able to handle many lanes with the limitations being
	system ram and hard drive space.
	
* The program only expects and recognizes 2 key files

* You must provide 2 key files even if they are identical

* intermediate files are written to the same output directory as summary data
	
* There is no check for the presence of an already modified second key file.  The 
	altered file needs to be manually removed after each call to the plugin
================================================================================
Version history

v1.1 = Added flags to modify read2 barcode, set stringeny on read1 and read2 barcodes,
	identify where barcodes and identifiers are used in the sequence files






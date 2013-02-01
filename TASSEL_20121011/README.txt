================================================================================
NCSU Fork of TASSEL for Genotyping by Sequencing

Software version: beta 0.9

Date: 2013-01-31

Legal information:
This software is provided "AS IS".  Use at your own discretion. 

Contact information:
Ross Whetten (ross_whetten@ncsu.edu)

================================================================================
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
* java 1.6+ (http://www.java.com)
* git 1.7+ (http://git-scm.com/)

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

** This is the bare bones approach to running only the paired end plugin. 
The TASSEL program can do many other things but requires a bit more setup
that can be read about at 
http://www.maizegenetics.net/tassel/docs/TasselPipelineGBS.pdf

--------------------------------------------------------------------------------
DOCUMENTATION:

* Fastq files must be named using the following template
	read#_flowcell_s_lane#_fastq.gz
	
	Working example:
	read1_C116JACXX_s_7_fastq.gz

* To run the paired end plugin, execute the following command (all on one line)
	via command line from the PROGRAM_DIR/.../outputs directory.  Optional arguements
	are in () and explained in GETTING STARTED TIPS.

	PROGRAM_DIR/.../maizegenetics/run_pairedEnd_pipeline.pl -[max_memory] 
	-fork1 -[plugin] -i [input_directory] -k [keyfile1:keyfile2] 
	-e [enzyme1:enzyme2] (-s [maximum_count]) (-c [minimum_detected]) 
	-o [output_directory] -endPlugin -runfork1

	Working example:
	PROGRAM_DIR/.../maizegenetics/run_pairedEnd_pipeline.pl -Xmx30g -fork1 
	-FastqPairedEndToTagCountPlugin -i fastq -k GBS1.key:GBS2.key 
	-e PstI-MspI:MspI-PstI -s 10000 -c 10 -o FastqPairedEndTagsAndTaxa 
	-endPlugin -runfork1
	
* The output files you WANT TO KEEP:
	* Barcode_ID_Info.txt = number of barcode combinations detected across all lanes
	* Summary_Stats.txt = printout of the number of forward, reverse, and pairs
		per lane over X interval
	* Tags_Totaled.txt = number of unique tags across all lanes
	* Tags_by_Taxa.txt = number of tags and taxa combination by lane over all lanes

* The plugin does create large intermediate files that can be deleted once the 
	program completes.  Any files in PROGRAM_DIR/outputs/FastqPariedEndTagsAndTaxa
	that start with "read...." and end with ".txt" can be removed to recover disk
	space.
	
* The plugin alters the second key file and stores it as [keyfile2]_mod.key.  
	Once processing of all data files is complete, this file can and should be 
	deleted.  If it is left within the directory it will be appended to each time
	the pluging is run.

--------------------------------------------------------------------------------
GETTING STARTED TIPS:

* Read through the entire README files before attempting to start

* Close any other programs to free system resources

* The program will take several hours to run, and will take longer with each lane
	of data that you process.
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
	
* Enzymes and keyfiles need to be separated by a ":".


================================================================================
Important known bugs, problems, or limitations:

2013-01-31: 
* Software has only been tested with 2 lanes of data (4 read files).
	In theory it should be able to handle many lanes with the limitations being
	system ram and hard drive space.
	
* The program only expects and recognizes 2 key files

* You must provide 2 key files even if they are identical

* intermediate files are written to the same output directory as summary data

* Summary_Stats.txt output doesn't label lane information.  Lanes are processed in
	numeric order so they are written to file that way, but it needs to be 
	explicitly added
	
* There is no check for the presence of an already modified second key file.  The 
	altered file needs to be manually removed after each call to the plugin
================================================================================
Version history







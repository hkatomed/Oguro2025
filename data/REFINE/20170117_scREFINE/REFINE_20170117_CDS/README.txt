===============================================================
README.txt	Information for Running REFINE Motif Prediction
===============================================================

1.	Software Requirements
2.	Setting up a REFINE Analysis
3.	Executing REFINE Script
4.	Output File Descriptions

===============================================================

1)	Software Requirements

	The following software programs/packages must be
	installed BEFORE running REFINE:

    A.  Perl (version 5.8.7 or later)is required.
	see http://www.cpan.org/
 
    B.  BioPerl (particularly the Seq and SeqIO modules)
	see http://www.bioperl.org/ 
    
    C.  R (version 2.1.1 or later)
	see http://www.r-project.org/

    D.  MEME (version 3.5.1)
	see http://meme.sdsc.edu/

    E.  DUST (Tatusov and Lipman, in preparation)
	see ftp.ncbi.nlm.nih.gov/pub/tatusov/dust

===============================================================

2)	Running REFINE

    A.  Be sure to install the required software described in Section 1.
    
    B.  The next step is to modify the file `refine.sh` so that it knows
	the locations of the required software programs from Section 1.
	Open refine.sh using a text editor program (such as emacs).
	You should see the following lines near the start of the file:

# IMPORTANT - USER MUST EDIT THESE PATHNAMES UPON INSTALLATION!!!
# Settings for Other Required Software
# 
set CMD_PERL = /usr/bin/perl            # modify this to point to Perl executable
set CMD_R    = /usr/local/bin/R         # modify this to point to R executable
set CMD_MEME = /usr/local/bin/meme      # modify this to point to MEME executable
set CMD_DUST = /usr/local/bin/dust.exe  # modify this to point to DUST executable

	Edit the pathname for each variable so that it points to the proper
	location of each executable program. (If you know you have installed R,
	for example, but are unsure of the path, you can type `which R` at the
	unix command line to find its location)
    
    C.  Next, setup the file(s) containing the reference sequences for motif analysis.
	If you are searching for motifs within S.cerevisiae in either the 200 base
	regions upstream (us) of the start codon or downstream (ds) of the stop codon
	for standard yeast genes, then you may use the provided non-redundant fasta-format
	sequence files yeast.nr.us.fa and yeast.nr.ds.fa in directory /REFINE_v0.1/seqs/.

	If you are using custom sequence files, make sure they are in fasta format and 
	contain sequences for all genes (both targets and non-targets).  For custom files,
	make sure the filenames end in either *.ds.fa or *.us.fa for compatibility, then
	place the custom sequence files in the directory /REFINE_v0.1./seqs/.
	
	For custom reference sequence files, you must also edit the following line in the
	file `refine.sh` to refer to the appropriate filename:

set REFSEQ = yeast.nr.$REG.fa        # filename for provided non-redundant yeast sequences
	
    	(i.e. change this to `set REFSEQ = custom.$REG.fa` if your new filename is custom.ds.fa)

    D.	Now, setup the file(s) containing background frequency values for mononucleotides
	present in the reference sequence database. If you are using the provided S.cerevisiae
	sequence files yeast.nr.us.fa or yeast.nr.ds.fa, then the accompanying provided files
	nt.freq.yeast.us or nt.freq.yeast.ds may be used as is. If you are using custom files,
	generate a text file in the format of nt.freq.yeast.ds that contains the background
	frequencies of each mononucleotide in the reference sequences. The file(s) should be
	placed in the directory /REFINE_v0.1/seqs/.

	Then edit the following line in `refine.sh` to refer to the custom file: 

set NTFREQ = nt.freq.yeast.$REG  # filename for provided background mononucleotide frequencies

    	(i.e. change this to `set REFSEQ = nt.custom.$REG` if your new filename is nt.custom.ds)

    E.  Prepare a text file listing the names of the target sequences (as named within
	the reference sequence file) with one sequence name per line, and save this with
	a filename such as "list.targets". (For using the pre-provided S. cerevisiae sequences,
	target genes should be specified using standard systematic yeast gene names, e.g. YFL039C).


===============================================================

    
3)	Executing REFINE

	Execute the REFINE analysis, for example by typing at the command line:
	
	`tcsh ./refine.sh targets ds`
	
	(corresponding to file list.targets and region 'ds')


===============================================================

4)	REFINE Output File Descriptions

	REFINE produces a number of text files as output. Assuming REFINE was
	executed with default options and the command `tcsh ./refine.sh targets ds`,
	the following files will be created in the subdirectory /output-targets/:
	
	A.	targets.ds.fa		
	
	This is a FASTA-format sequence file containing the 
	subset of sequences present in the reference sequence file 
	that are specified as targets in the file "list.targets".
	
	B.	rel.sp.6.targets.ds.txt	
	
	This is a text file with information about the relative enrichment of
	nucleotide hexamers in the target sequences versus the reference set,
	with each line in the following format:

	AAAAAA	14	33	2294	5675	-0.325334912174749
	AAAAAC	9	33	1202	5675	-0.598896497462366
	AAAAAG	10	33	1916	5675	-0.140795589963578
	AAAAAT	10	33	2074	5675	-0.0847664575737535

	column-1:	sequence of nucleotide hexamer
	column-2:	number of target sequences containing this hexamer
	column-3:	number of total target sequences evaluated
	column-4:	number of reference sequences containing this hexamer
	column-5:	number of total reference sequences evaluated
	column-6:	-log10(p-value) for significance of observed relative enrichment,
			based on hyper-geometric distribution

	C.	wds.3.6.targets.ds.txt	
	
	This is a text file containing a list of hexamers with significant enrichment,
	in the targets relative to the rest of the reference set (default p<=0.001).
	
	D.	mask.3.targets.ds.fa
	
	This is a FASTA-format file containing filtered target sequences, which is
	used as input for MEME motif analysis. 
		
	E.	meme.mask.3.targets.ds
	
	This MEME-generated output file has information on up to three predicted motifs.
	
	
	The remaining output file types (F-I) each correspond to a single predicted motif,
	with 'X' indicating the numbering of motifs (X = 1 or 2 or 3):
	
	F.	blk.targets.ds.X
	
	This text file lists sequences for each instance of motif X identified by MEME.
		
	G.	mls.targets.ds.X.txt
	
	This is a text file with information on the highest-scoring possible motif site
	for every sequence in the reference set in the following format:

	YAL001C	1.071092	24	34	TTTTATTTCTT	0
	YAL002W	0.701027	181	191	CATTACTCATT	0
	YAL003W	1.581643	145	155	CCTCATGCTTT	0
	YAL005C	0.803232	90	100	TCCCAATTCTT	0
	
	column-1:	Name of the sequence from the reference file
	column-2:	Motif score (based on parameters in file lod.targets.ds.X)
	column-3:	Starting coordinate of the site
	column-4:	Ending coordinate of the site
	column-5:	Sequence of the site
	column-6:	1 for target sequences from list.targets, 0 otherwise
		
	H.	lod.targets.ds.X
	
	This file contains the log-odds scoring matrix data used for evaluating
	and classifying potential motif sites, in the following format:
	
	CUTOFF= 3.754213
	ALPHABET= ACGT
	log-odds matrix: alength= 4 w= 11
	0.302264	-0.065306	-1.017758	-0.317144	
	-0.384306	-0.348607	0.184662	0.219099	
	-0.639579	0.274642	-0.346816	0.219099	
	-1.338549	0.708298	-1.045786	-0.600445	
	0.446781	-1.047577	-1.045786	-1.299415	
	-0.639579	-0.348607	-1.045786	0.424861	
	-1.338549	-1.047577	-1.045786	0.485915	
	-1.338549	0.737753	-1.045786	-1.299415	
	-0.639579	0.642619	-1.045786	-0.345172	
	-0.384306	-1.047577	-1.045786	0.424861	
	-0.384306	-1.047577	-1.045786	0.424861	
	
	The 'CUTOFF' value indicates the score threshold required for
	positively classifying a potential site as a motif.
	
	The 'w' value indicates the length of the motif (# of matrix rows).
	
	The four numerical values in each row (1 through w) indicate
	the value associated with observing an A,C,G, or U at the
	corresponding position (1 through w) in the potential motif site.
	The first row of the matrix corresponds to the first base in the
	motif, and so on...
	
	I.	sites.targets.ds.X.txt	
	
	This text file contains information on sites in the reference
	sequence database that are classified as motif instances
	(i.e. the score for the site based on the log-odds scoring matrix
	exceeds the cutoff value for motif classification).
	
	Columns follow the same format as for mls.targets.ds.X.txt:

	column-1:	Name of the sequence from the reference file
	column-2:	Motif score (based on parameters in file lod.targets.ds.X)
	column-3:	Starting coordinate of the site
	column-4:	Ending coordinate of the site
	column-5:	Sequence of the site
	column-6:	1 for target sequences from list.targets, 0 otherwise


===============================================================
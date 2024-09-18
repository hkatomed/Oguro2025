#!/bin/tcsh
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# Shell Script for Running REFINE Motif Prediction
#
# Please see README.txt for information about using this script!
#

# IMPORTANT - USER MUST EDIT THESE VARIABLES UPON INSTALLATION!!!
# Settings for Other Required Software
# 
set CMD_PERL = /usr/bin/perl            # modify this to point to Perl executable
set CMD_R    = /usr/bin/r         # modify this to point to R executable
set CMD_MEME = /Users/hrk_kato/meme/bin/meme      # modify this to point to MEME executable
set CMD_DUST = /Users/hrk_kato/meme/bin/dust  # modify this to point to DUST executable
set CMD_FASTAGETMARKOV = /Users/hrk_kato/meme/bin/fasta-get-markov


# Input Arguments (i.e. `tcsh refine.sh rbp-name ds`)
set NAME   = $1                  # input argument is suffix of file in current directory list.NAME
set REG    = $2                  # region for sequence analysis (ds or us)
set PREFIX = $3
set REFSEQ = $PREFIX.$REG.fa    # filename for provided non-redundant yeast sequences
set NTFREQ = nt.freq.$PREFIX.$REG  # filename for provided background mononucleotide frequencies


# Directory Settings
set DIR_HOME = `pwd`
set DIR_OUT  = $DIR_HOME/output-$NAME-$REG-$PREFIX-K$4-$5
set DIR_MEMEOUT  = $DIR_HOME/meme_out-$NAME-$REG-$PREFIX-K$4-$5
set DIR_SEQ  = $DIR_HOME/seqs  # directory for sequence files used in REFINE
set DIR_PERL = $DIR_HOME/src  # directory for additional scripts used in REFINE
set DIR_TMP  = $DIR_HOME/temp  # directory for temporary files


# Parameters for REFINE Prediction
set K = $4      # nucleotide length (use 6 for hexamers)
set C = 3      # -log10(p-value) cutoff for significant hexamers (i.e. p<= 0.001)
set FLANK = 3  # number of bases adjacent to significant k-mers to be retained
set GAP = 6    # maximum distance of gaps between significant k-mers to be retained
set CT = 3     # minimum number of significant k-mer sites for target sequences to be kept
set XMAX = 15  # maximum length of consecutive X-bases (for collapsing filtered segments)

# Parameters for MEME
set MINW = 6          # minimum motif width
set MAXW = 15         # maximum motif width
set NMOTIFS = 3       # max number of motifs to be reported
set MAXSIZE = 20000000  # max size of input sequence
set EVT  = 10         # E-value threshold


echo "Running REFINE for target genes in list.$NAME for $REG region"
if (! -d $DIR_OUT) then
    mkdir $DIR_OUT
endif

if (! -d $DIR_TMP) then
    mkdir $DIR_TMP
endif

# Preparation of Target Sequence Files
$CMD_PERL $DIR_PERL/selectfaseqs.pl -list $DIR_HOME/list.$NAME -i $DIR_SEQ/$REFSEQ -o $DIR_OUT/tmp.$NAME.$REG.fa
echo "selectfaseqs.pl DONE"
$CMD_PERL $DIR_PERL/fa2longfa.pl -i $DIR_OUT/tmp.$NAME.$REG.fa > $DIR_OUT/$NAME.$REG.fa
echo "fa2longfa.pl DONE"
$CMD_DUST $DIR_OUT/$NAME.$REG.fa > $DIR_OUT/tmp.dust.$NAME.$REG.fa
echo "dust DONE"
$CMD_PERL $DIR_PERL/fa2longfa.pl -i $DIR_OUT/tmp.dust.$NAME.$REG.fa > $DIR_OUT/dust.$NAME.$REG.fa
echo "fa2longfa.pl DONE"

# Identification of Significantly Enriched Hexamers
$CMD_PERL $DIR_PERL/seqPres.pl -ss -size $K $DIR_OUT/$NAME.$REG.fa > $DIR_OUT/sp.$K.$NAME.$REG
$CMD_PERL $DIR_PERL/seqPres.pl -ss -size $K $DIR_SEQ/$REFSEQ > $DIR_SEQ/sp.$K.all.$REG
$CMD_PERL $DIR_PERL/splice.pl $DIR_OUT/sp.$K.$NAME.$REG $DIR_SEQ/sp.$K.all.$REG | cut -f1-3,5-6 > $DIR_TMP/tempfile.txt
cp $DIR_PERL/*.R $DIR_TMP/.
cd $DIR_TMP
$CMD_R CMD BATCH subsethyperpvals.R
mv $DIR_TMP/tempfile.p $DIR_OUT/rel.sp.$K.$NAME.$REG.txt
$CMD_PERL $DIR_PERL/cutoff.pl -f6 -lo -$C $DIR_OUT/rel.sp.$K.$NAME.$REG.txt | cut -f1 | sort -u > $DIR_OUT/wds.$C.$K.$NAME.$REG.txt
rm $DIR_TMP/tempfile.txt
cd $DIR_HOME

# Masking of Selected Target Sequences
$CMD_PERL $DIR_PERL/wordmask.pl -ss -x -flank $FLANK -gap $GAP -ct $CT -max $XMAX -wds $DIR_OUT/wds.3.$K.$NAME.$REG.txt -i $DIR_OUT/dust.$NAME.$REG.fa -o $DIR_OUT/mask.$C.$NAME.$REG.fa

# MEME Prediction of Motifs within Masked Target Sequences
if($5 == 0th) then
echo "Markov: 0th"
$CMD_MEME -oc $DIR_MEMEOUT $DIR_OUT/mask.$C.$NAME.$REG.fa -dna -minw $MINW -maxw $MAXW -maxsize $MAXSIZE -evt $EVT -nmotifs $NMOTIFS > $DIR_OUT/meme.mask.$C.$NAME.$REG
$CMD_MEME -oc $DIR_MEMEOUT $DIR_OUT/mask.$C.$NAME.$REG.fa -dna -minw $MINW -maxw $MAXW -text -maxsize $MAXSIZE -evt $EVT -nmotifs $NMOTIFS > $DIR_OUT/meme.mask.$C.$NAME.$REG

else if($5 == 5th) then
echo "Markov: 5th"
echo "Checking preexisting markov_5.txt"
if(-f $DIR_SEQ/markov_5.txt) then
rm $DIR_SEQ/markov_5.txt
echo "Preexisting markov_5.txt deleted"
else
echo "markov_5.txt not exist."
endif

echo "Creating markov_5.txt"
$CMD_FASTAGETMARKOV -m 5 $DIR_SEQ/$REFSEQ > $DIR_SEQ/markov_5.txt
if(-f $DIR_SEQ/markov_5.txt) then
echo "markov_5.txt created."
else
echo "Error: markov_5.txt not found."
echo ""
endif
$CMD_MEME -oc $DIR_MEMEOUT $DIR_OUT/mask.$C.$NAME.$REG.fa -dna -minw $MINW -maxw $MAXW -maxsize $MAXSIZE -evt $EVT -nmotifs $NMOTIFS -bfile $DIR_SEQ/markov_5.txt > $DIR_OUT/meme.mask.$C.$NAME.$REG
$CMD_MEME -oc $DIR_MEMEOUT $DIR_OUT/mask.$C.$NAME.$REG.fa -dna -minw $MINW -maxw $MAXW -text -maxsize $MAXSIZE -evt $EVT -nmotifs $NMOTIFS -bfile $DIR_SEQ/markov_5.txt > $DIR_OUT/meme.mask.$C.$NAME.$REG
else
echo "Error: The 5th argument (Markov) should be 0th or 5th."
endif

# Process Individual Predicted Motifs
foreach k (1 2 3)
$CMD_PERL $DIR_PERL/meme2blocks.pl -n $k $DIR_OUT/meme.mask.$C.$NAME.$REG > $DIR_OUT/blk.$NAME.$REG.$k

if (-z $DIR_OUT/blk.$NAME.$REG.$k) then
rm $DIR_OUT/blk.$NAME.$REG.$k
else

$CMD_PERL $DIR_PERL/getlodfile.pl -i $DIR_OUT/blk.$NAME.$REG.$k -bfile $DIR_SEQ/$NTFREQ -base 10 -pseud 0.25 > $DIR_OUT/tmp.lod.$NAME.$REG.$k
$CMD_PERL $DIR_PERL/getmaxlodsite.pl -lod $DIR_OUT/tmp.lod.$NAME.$REG.$k -i $DIR_SEQ/$REFSEQ -o $DIR_OUT/mls.$NAME.$REG.$k.txt
$CMD_PERL $DIR_PERL/membership.pl -list $DIR_HOME/list.$NAME -i $DIR_OUT/mls.$NAME.$REG.$k.txt > $DIR_TMP/tempfile.txt
cd $DIR_TMP
$CMD_R CMD BATCH scoresites.R
cd $DIR_HOME

mv $DIR_TMP/tempfile.ss $DIR_OUT/tmp.ss.$NAME.$REG.$k
mv $DIR_TMP/tempfile.txt $DIR_OUT/mls.$NAME.$REG.$k.txt
$CMD_PERL $DIR_PERL/updatelodfile.pl $DIR_OUT/tmp.ss.$NAME.$REG.$k $DIR_OUT/tmp.lod.$NAME.$REG.$k > $DIR_OUT/lod.$NAME.$REG.$k
$CMD_PERL $DIR_PERL/getalllodsites.pl -lod $DIR_OUT/lod.$NAME.$REG.$k -i $DIR_SEQ/$REFSEQ -o $DIR_OUT/sites.$NAME.$REG.$k.txt
$CMD_PERL $DIR_PERL/membership.pl -list $DIR_HOME/list.$NAME -i $DIR_OUT/sites.$NAME.$REG.$k.txt > $DIR_OUT/temp
mv $DIR_OUT/temp $DIR_OUT/sites.$NAME.$REG.$k.txt
endif
end

# Removal of temporary files
rm $DIR_OUT/tmp.*
rm $DIR_OUT/sp.$K.$NAME.$REG
rm $DIR_OUT/dust.$NAME.$REG.fa
rm -rf $DIR_TMP/

echo ""
echo "REFINE has completed."
echo "Output files are in $DIR_OUT/"
echo ""

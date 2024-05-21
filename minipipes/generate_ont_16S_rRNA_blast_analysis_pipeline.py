#!/usr/bin/perl -w

# The MIT License
# Copyright (c) 2021 Adrian Tan <adrian_tan@nparks.gov.sg>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

generate_nanopore_16S_rRNA_analyses_pipeline

=head1 SYNOPSIS

 generate_nanopore_16S_rRNA_blast_analysis_pipeline [options]
  -s     sample file
         column 1: sample ID
         column 2: nanopore barcode
  -o     output directory
  -m     make file name

=head1 DESCRIPTION

This script implements the pipeline for mapping the reads against a 16S rRNA database.

=cut

my $help;
my $sampleFile;
my $outputDir;
my $makeFile = "nanopore_16S_rRNA_blast_analysis.mk";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help,
                's=s'=>\$sampleFile,
                'o=s'=>\$outputDir,
                'm:s'=>\$makeFile
               )
  || !defined($outputDir)
  || !defined($sampleFile))
{
    if ($help)
    {
        pod2usage(-verbose => 2);
    }
    else
    {
        pod2usage(1);
    }
}

$makeFile = "$outputDir/$makeFile";

#programs
my $blastn= "/usr/local/ncbi-blast-2.11.0+/bin/blastn";
my $seqtk = "/usr/local/seqtk-1.3/seqtk";
my $nanoplot = "/usr/local/bin/NanoPlot";

printf("generate_16S_rRNA_blast_analysis_pipeline.pl\n");
printf("\n");
printf("options: output dir           %s\n", $outputDir);
printf("         make file            %s\n", $makeFile);
printf("         sample file          %s\n", $sampleFile);
printf("\n");

################################################
#Helper data structures for generating make file
################################################
my @tgts = ();
my @deps = ();
my @cmds = ();
my $tgt;
my $dep;
my $log;
my $err;
my @cmd;

#################
#Read sample file
#################
my %SAMPLE = ();
my @SAMPLE = ();

open(SA,"$sampleFile") || die "Cannot open $sampleFile\n";
while (<SA>)
{
    s/\r?\n?$//;
    if(!/^#/)
    {
        my ($sampleID, $fastqFile) = split(/\s+/, $_);

        if (exists($SAMPLE{$sampleID}))
        {
            exit("$sampleID already exists. Please fix.");
        }

        $SAMPLE{$sampleID} = ();
        $SAMPLE{$sampleID}{FASTQ} = $fastqFile;
        push(@SAMPLE, $sampleID);
    }
}
close(SA);
print "read in " . scalar(@SAMPLE) . " samples\n";

#####################################
#convert to fasta, filter then sample
#####################################
for my $sampleID (@SAMPLE)
{
    my $sampleDir = "$outputDir/samples/$sampleID";
    mkpath($sampleDir);
    my $outputFASTAFile = "$sampleDir/filtered.sampled.fa";

    $tgt = "$outputFASTAFile.OK";
    $dep = "";
    $log = "$outputFASTAFile.log";
    $err = "$outputFASTAFile.err";
    @cmd = ("$seqtk seq -a -L 1400 $SAMPLE{$sampleID}{FASTQ} | $seqtk sample - -s 13 10000 > $outputFASTAFile");
    makeJob("local", $tgt, $dep, @cmd);
}

############################
#construct 16S rRNA database
############################
#esearch -db nucleotide -query "33175[BioProject] OR 33317[BioProject] " | efetch -format fasta > out.fasta
#makeblastdb -in 21940seq_16s.fasta  -dbtype nucl -parse_seqids
#export BLASTDB=/usr/local/ncbi-blast-2.11.0+/bin

################################
#blast against 16S rRNA database
################################
for my $sampleID (@SAMPLE)
{
    my $sampleDir = "$outputDir/samples/$sampleID";
    my $inputFASTAFile = "$sampleDir/filtered.sampled.fa";
    my $outputBlastFile = "$sampleDir/blast_results.txt";

    $tgt = "$outputBlastFile.OK";
    $dep = "$inputFASTAFile.OK";
    $log = "$outputBlastFile.log";
    $err = "$outputBlastFile.err";
    @cmd = ("blastn -db $outputDir/db/21940seq_16s.fasta -query $inputFASTAFile -outfmt \"6 stitle pident\" -out $outputBlastFile -max_target_seqs 1 2> $err > $log");
    makeJob("local", $tgt, $dep, @cmd);
}

########################
#summarise blast results
########################
for my $sampleID (@SAMPLE)
{
    my $sampleDir = "$outputDir/samples/$sampleID";
    my $inputTXTFile = "$sampleDir/blast_results.txt";
    my $outputTXTFile = "$sampleDir/blast_results.summarised.txt";

    $tgt = "$outputTXTFile.OK";
    $dep = "$inputTXTFile.OK";
    $log = "$outputTXTFile.log";
    $err = "$outputTXTFile.err";
    @cmd = ("cat $inputTXTFile | cut -f1  | perl -lane '{/([^\s]+) ([^\s]+)/; print \"\$1 \$2\"}' | sort | uniq -c | sort -nrk1 > $outputTXTFile");
    makeJob("local", $tgt, $dep, @cmd);
}

#*******************
#Write out make file
#*******************
GENERATE_MAKEFILE:
print "\nwriting makefile\n";

open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
print MAK ".DELETE_ON_ERROR:\n\n";
print MAK "all: @tgts\n\n";

#clean
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $outputDir/*.* $outputDir/intervals/*.*");

for(my $i=0; $i < @tgts; ++$i) {
    print MAK "$tgts[$i] : $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#Functions
##########

#run a job either locally or by pbs
sub makeJob
{
    my ($method, @others) = @_;

    if ($method eq "local")
    {
        my ($tgt, $dep, @rest) = @others;
        makeLocalStep($tgt, $dep, @rest);
    }
    else
    {
        die "unrecognized method of job creation : $method\n";
    }
}

sub makeLocalStep
{
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        if ($cmd =~ /\|/)
        {
            $cmd .= "\tset -o pipefail; " . $c . "\n";
        }
        else
        {
            $cmd .= "\t$c\n";
        }
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}

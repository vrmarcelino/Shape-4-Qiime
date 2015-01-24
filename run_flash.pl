#!/usr/bin/perl
@ff = <*_R1_001.fastq_trimmed.fastq>;

foreach $file1 (@ff) {
	$file1 =~ /^(.*)_R1_001.fastq_trimmed.fastq/;
	$root = $1;
	$file2 = $root."_R2_001.fastq_trimmed.fastq";
	$outfile = $root."_paired.fastq";
	#print "$file1 + $file2 => $outfile\n";
	system("flash $file1 $file2 -m 70 -M 305 -O -x 0.1 -o $outfile -d Merged")
}


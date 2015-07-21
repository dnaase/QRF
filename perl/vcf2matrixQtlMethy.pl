 # The MIT License (MIT)
 # Copyright (c) 2015 dnaase <Yaping Liu: lyping1986@gmail.com>

 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:

 # The above copyright notice and this permission notice shall be included in all
 # copies or substantial portions of the Software.

 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.

#!/usr/bin/perl -w

## vcf2matrixQtlFormat.pl
##convert Bis-SNP output cpg.vcf into a DNA methylaiton input for matrix QTL software.

## Author: yaping  lyping1986@gmail.com
## Nov 26, 2014:10:44:15 AM

use strict;
use Getopt::Long;
use File::Basename;

sub usage{
	print STDERR "\nUsage:\n";
    print STDERR "perl vcf2matrixQtlMethy.pl [Options] input_file_name [CG] \n\n";
    print STDERR "convert Bis-SNP output cpg.vcf into a DNA methylaiton input for matrix QTL software.\n\n";

	
	print STDERR "  [Options]:\n\n";	
 	print STDERR "  --start_sample NUM: the start sample number in VCF file to output. (Default: not enable. Output all of samples)\n\n";
 	print STDERR "  --end_sample NUM: the end sample number in VCF file to output. (Default: not enable. Output all of samples)\n\n";
 	exit(1);
}

##default option setting
my $start_sample = "";
my $end_sample = "";

GetOptions( "start_sample=i" => \$start_sample,
			"end_sample=i" => \$end_sample,
);

usage() if ( scalar(@ARGV) == 0 );

my $input_file_name = $ARGV[0];
my $type = $ARGV[1] || "\\w+";

my $cpg_name_output = $input_file_name;
$cpg_name_output =~ s/\.\w+$/.matrixQtlInput.txt/;
open(OUT,">$cpg_name_output") or die "can't write file $cpg_name_output:$!\n";


#if($type eq ""){
#	$type = "\\w+"; ##by default, output all sites in VCF file
#}
my $sample_num=0;
my $cpg_num=0;
open(FH,"<$input_file_name") or die "can't open file $input_file_name :$!\n";
while(<FH>){
	my $line=$_;
	chomp($line);
	next if $line =~ /^##/;
	my @splitin = split "\t", $line;
	
	##header
	if($splitin[0] eq "#CHROM"){
		print OUT "id";
		##output sample number
		if($start_sample eq ""){
			$start_sample=9;
		}else{
			$start_sample+=8;
		}
		if($end_sample eq ""){
			$end_sample=$#splitin;
		}else{
			$end_sample+=8;
		}
		for(my $i=$start_sample;$i<=$end_sample;$i++){
			print OUT "\t$splitin[$i]";
			$sample_num++;	
		}
		print OUT "\n";
	}
	
	##main body
	next unless ($splitin[6] eq "PASS");
	print OUT "cg:$splitin[0]-$splitin[1]";
	for(my $i=$start_sample;$i<=$end_sample;$i++){
			#print "$start_sample\t$end_sample\t$i\n";
			if($splitin[$i] eq "./."){
				print OUT "\tNA";
			}elsif($splitin[$i] =~ /:(\d+):$type:(\d+):/){
				my $num_c = $1;
				my $num_t = $2;
				if($num_c + $num_t > 0){
					my $methy = $num_c/($num_c + $num_t);
					$methy = sprintf("%.5f",$methy);
					print OUT "\t$methy";
				}else{
					print OUT "\tNA";
				}
			}else{
				print OUT "\tNA";
			}
	}
	print OUT "\n";
	$cpg_num++;
}

close(FH);
close(OUT);
print STDERR "Output $sample_num samples and $cpg_num lines in total\n\n";

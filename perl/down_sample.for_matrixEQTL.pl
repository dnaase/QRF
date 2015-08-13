#!/usr/bin/perl -w

use strict;

my $configure_file=$ARGV[0] or die "no input configure files!\n";
my $sample_size=$ARGV[1] or die "no down sample size!\n";
my $configure_file_down_sampled=$ARGV[2] or die "no output downsampled configure file name!\n";

my $qrf_path=`echo \$QRF`;
chomp($qrf_path);

#get configure parameters
my ($dir,$snp_tab,$gene_tab,$cov_tab) = check_conf_file($configure_file);

#downsample and get new SNP info, gene info and cov,info
my $cmd="R --no-save --no-restore --args wd=$dir snpInfo=$snp_tab exprInfo=$gene_tab covarInfo=$cov_tab sampleSize=$sample_size < $qrf_path/R/downsample.R\n";
run_cmd($cmd);

#generate new configure files
open(FH,"<$configure_file") or die "can't open configure file $configure_file:$!\n";
open(OUT,">$configure_file_down_sampled") or die "can't open new configure file $configure_file_down_sampled:$!\n";
while(<FH>){
		chomp;
		my @f=split "=";
		
		if($f[0] eq "SNP_tab"){
			$snp_tab=$f[1].".SampleSize-${sample_size}.txt";
		}elsif($f[0] eq "gene_tab"){
			$gene_tab=$f[1].".SampleSize-${sample_size}.txt";
		}elsif(($f[0] eq "covar_tab") and ($f[1] ne "NULL") and ($f[1] ne "null") and ($f[1] ne "")){
			$cov_tab=$f[1].".SampleSize-${sample_size}.txt";
		}
		print OUT join("=",@f)."\n";
}
close(FH);
close(OUT);



sub check_conf_file{
	my $input=shift @_;
	my $dir="";
	my $snp_tab="";
	my $gene_tab="";
	my $cov_tab="";
	
	open(FH,"<$input") or die "can't open configure file $input:$!\n";
	while(<FH>){
		chomp;
		my @f=split "=";
		if($f[0] eq "directory"){
			$dir=$f[1];
		}elsif($f[0] eq "SNP_tab"){
			$snp_tab=$f[1];
		}elsif($f[0] eq "gene_tab"){
			$gene_tab=$f[1];
		}elsif($f[0] eq "covar_tab"){
			$cov_tab=$f[1];
		}
	}
	close(FH);
	
	##check parameter read
	if(($dir eq "") or ($dir eq "NULL")){
		$dir="./";
	}
	if(($snp_tab eq "") or ($snp_tab eq "NULL") or ($gene_tab eq "") or ($gene_tab eq "NULL")){
		die "SNP information and Gene/CpG information could not be NULL!!\n\n";
	}
	
	return($dir,$snp_tab,$gene_tab,$cov_tab);
}


sub run_cmd{
	my $cmd=shift @_;
	print STDERR "$cmd\n";
	system($cmd) == 0 || die "can't execute command $cmd:$!\n";
}



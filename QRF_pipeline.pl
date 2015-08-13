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
##QRF_pipeline.pl: the main entry of QRF software


## author: Yaping Liu  lyping1986@gmail.com 
## time: 2015-7-21

#Usage: perl QRF_pipeline.pl [option] configure.txt

use strict;
use Getopt::Long;
use File::Basename;


sub usage {

    print STDERR "\nUsage:\n";
    print STDERR "perl QRF_pipeline.pl [option] prefix configure.txt\n\n";
    print STDERR "Require R-3.0 (MatrixEQTL library) or above, Java-1.6 or above and bigWigAverageOverBed (download from UCSC utils) to installed\n\n";

    print STDERR " [Options]:\n\n";
    print STDERR " ###############General options:\n\n";
	print STDERR "  --mode NUM : 1. Detect meQTL/eQTL;(Default mode)\n";
	print STDERR "  			 2. Generate training file;\n";
	print STDERR "  --region STR : Specify the region for mode 1 or 2. Should be format like: chr1:1-1000 or chr1 (Default: chr1).\n";
	print STDERR "  --fdr NUM : Specify the false discovery rate for Matrix EQTL and QRF (Default: 0.01).\n";
	print STDERR "  --mem NUM : Specify the number of gigabytes in the memory to use (Default: 15).\n";
	print STDERR "  --qrf_path DIR : Specify the QRF root direcotry (Default: not specified. use environment variable \$QRF).\n\n";

    print STDERR " ###############Options for mode 1:\n\n";
	print STDERR "  --tree_num NUM : Specify the number of tree used for QRF (Default: 1000).\n";
	print STDERR "  --class_index NUM : Specify the column number that are the label column for QRF (Default: 5).\n";
	print STDERR "  --label_class STR : Specify the class label name used for QRF (Default: meqtl).\n";
	print STDERR "  --permutation_times NUM : Specify the number of permutation used for QRF. 0 or negative value means not enabled (Default: 1).\n";
	print STDERR "  --sep STR : Specify the string used to seperate the column (Default: \\\\t, tab delimit).\n";
	print STDERR "  --hic_resolution NUM : Specify the resolution of HiC signal used, default is 1kb (Default: 1).\n\n";

    print STDERR " ###############Options for mode 2:\n\n";
	print STDERR "  --positive_probes_sampling NUM : Specify the number of positive probes used for QRF training (Default: 10000).\n";
	print STDERR "  --negative_probes_sampling NUM : Specify the number of negative probes used for QRF training (Default: 10000).\n\n";

    print STDERR " ###############Options for mode 3:\n\n";
	print STDERR "  --r2 NUM : Specify the r square threshold for LD independent interval generation (Default: 0.8).\n";
	print STDERR "  --ld_dist NUM : Specify the distance threshold for LD independent interval generation (Default: 250000).\n\n";

			
    exit(1);
}


##default option setting
my $mode=1;

my $mem=15;
my $region="chr1";
my $fdr=0.01;
my $qrf_path=`echo \$QRF`;
chomp($qrf_path);

my $tree_num=1000;
my $class_index=4;
my $label_class="meqtl";
my $permutation_times=1;
my $sep="\\\\t";
my $hic_resolution=1;

my $positive_probes_sampling=10000;
my $negative_probes_sampling=10000;

my $omit_matrixeqtl="";
my $omit_matrixeqtl_result_anno="";
my $omit_add_hic="";
my $omit_add_recomb="";
my $omit_cal_recomb_oe="";
my $omit_make_input_matrix="";
my $omit_random_forest="";
my $omit_qrf_result_anno="";

GetOptions( 
			"mode=i" => \$mode,
			"mem=i" => \$mem,
			"region=s" => \$region,
			"fdr=f" => \$fdr,
			"qrf_path=s" => \$qrf_path,
			"tree_num=i" => \$tree_num,
			"class_index=i" => \$class_index,
			"label_class=s" => \$label_class,
			"permutation_times=i" => \$permutation_times,
			"sep=s" => \$sep,
			"hic_resolution=i" => \$hic_resolution,
			"positive_probes_sampling=i" => \$positive_probes_sampling,
			"negative_probes_sampling=i" => \$negative_probes_sampling,

			"omit_matrixeqtl" => \$omit_matrixeqtl,
			"omit_matrixeqtl_result_anno" => \$omit_matrixeqtl_result_anno,
			"omit_add_hic" => \$omit_add_hic,
			"omit_add_recomb" => \$omit_add_recomb,
			"omit_cal_recomb_oe" => \$omit_cal_recomb_oe,
			"omit_make_input_matrix" => \$omit_make_input_matrix,
			"omit_random_forest" => \$omit_random_forest,
			"omit_qrf_result_anno" => \$omit_qrf_result_anno,

);



check_parameter();

my $prefix=$ARGV[0];
my $conf_file=$ARGV[1];

my ($dir,$matrixqtlresult,$snp_loc,$gene_loc,$snp_tab,$gene_tab,$cov_tab,$recomb, $hic, $train, $recombination_expect) = check_conf_file($conf_file);

my $chr=$region;
if($chr=~/\b\w+:\S+/){
		$chr=~s/:\S+//;
		
}

##call meQTL/eQTL by MatrixEQTL
if(($omit_matrixeqtl eq "") and (($matrixqtlresult eq "") or ($matrixqtlresult eq "NULL"))){

	my $tmp=`wc $dir/$snp_tab`;
	chomp($tmp);
	my @f=split " ",$tmp;
	my $sampleSize=int($f[1]/$f[0])-1;
	print STDERR "no Matrix EQTL results provided. So called MatrixEQTL now ...\n\n";
	my $cmd="R --no-save --no-restore --args qrf_path=$qrf_path wd=$dir snpInfo=$snp_tab exprInfo=$gene_tab  snpLoc=$snp_loc exprLoc=$gene_loc covarInfo=$cov_tab prefix=$prefix chr=$chr sampleSize=$sampleSize < $qrf_path/R/call_meqtl_random.byChr.R\n";
	run_cmd($cmd);
	$matrixqtlresult = "$dir/cis-qtl.matrixEQtlAll.SampleSize-$sampleSize.$prefix.$chr.txt";
}else{
	$matrixqtlresult="$dir/$matrixqtlresult";
}

##annotate MatrixEQTL result with coordinate
my $matrixqtlresult_cor=$matrixqtlresult;
$matrixqtlresult_cor=~s/\.\w+$/.sig.sameChr.addCor.bed/;
my $matrixqtlresult_cor_uniq=$matrixqtlresult_cor;
$matrixqtlresult_cor_uniq=~s/\.\w+$/.uniq.bed/;
my $cmd = "perl $qrf_path/perl/add_cor_to_matrixqtl.pl $matrixqtlresult $dir/$snp_loc $dir/$gene_loc\n";
if($omit_matrixeqtl_result_anno eq ""){
	print STDERR "Adding genomic coordinate to MatrixEQTL result ...\n\n";
	run_cmd($cmd);
	#get uniq bed file
	$cmd = "perl $qrf_path/perl/uniqLine.pl $matrixqtlresult_cor $matrixqtlresult_cor_uniq\n";
	run_cmd($cmd);
	
}



##attach HiC signal
my $matrixqtlresult_hic=$matrixqtlresult_cor_uniq;
$matrixqtlresult_hic=~s/\.\w+$/.hic_signal.bed/;
if($omit_add_hic eq ""){
	print STDERR "Adding HiC signal to MatrixEQTL result ...\n\n";
	$cmd = "perl $qrf_path/perl/get_hic_freq.sparse.pl $hic_resolution $matrixqtlresult_cor_uniq $dir/$hic 1 1 > $matrixqtlresult_hic\n";
	run_cmd($cmd);
	
}

##attach recombination rate, and calculate O-E

my $matrixqtlresult_recomb_prefix=$matrixqtlresult_cor_uniq;
$matrixqtlresult_recomb_prefix=~s/\.\w+$/.recomb_signal/;
if($omit_add_recomb eq ""){
	print STDERR "Adding recombination rate signal to MatrixEQTL result ...\n\n";
	$cmd = "perl $qrf_path/perl/alignWig2Bed.bigWigAverageOverBed.pl $matrixqtlresult_recomb_prefix $matrixqtlresult_cor_uniq $dir/$recomb --min_data 1 \n";
	run_cmd($cmd);
	
}
my $recomb_prefix=basename($recomb);
$recomb_prefix=~s/(\w+)\S+/$1/;
my $matrixqtlresult_recomb=$matrixqtlresult_recomb_prefix.".alignTo.$recomb_prefix.min_data-1.txt";

my $matrixqtlresult_recomb_oe=$matrixqtlresult_recomb;
$matrixqtlresult_recomb_oe=~s/\.\w+$/.recomb_observed_expect.txt/;
if($omit_cal_recomb_oe eq ""){
	print STDERR "Calculate Observed-Expectation recombination rate signal to MatrixEQTL result ...\n\n";
	$cmd = "perl $qrf_path/perl/calculate_recomb_OE.from_summary_random.pl $matrixqtlresult_recomb $dir/$recombination_expect > $matrixqtlresult_recomb_oe\n";
	run_cmd($cmd);
	
}

##make matrix file for random forest training
my $input_matrix=$matrixqtlresult_cor_uniq;
$input_matrix=~s/\.\w+$/.tstat_hic_recomb_oe.txt/;
if($omit_make_input_matrix eq ""){
	$cmd="paste $matrixqtlresult_hic $matrixqtlresult_recomb_oe | cut -f7,9-10,14 | perl -ne 'chomp;\@f=split \"\\t\";\$f[0]=abs(\$f[0]);if(\$f[1]<$fdr){\$name=\"$label_class\";}else{\$name=\"no$label_class\";}print \"\$f[0]\\t\$f[2]\\t\$f[3]\\t\$name\\n\";' > $input_matrix\n";
	run_cmd($cmd);
	
}

##run random forest model
my $output_matrix=$input_matrix;
$output_matrix=~s/\.\w+$/.afterRandomForest.txt/;
if($omit_random_forest eq ""){
	print STDERR "QRF core model ...\n\n";
	$cmd = "java -Xmx${mem}G -jar $qrf_path/dist/QRF.jar $dir/$train -sep $sep -classIndex $class_index -outputFile $output_matrix -inputFile $input_matrix -ignoreCV -numTrees $tree_num -permutation -permutationNum $permutation_times -permutationClass $label_class\n";
	run_cmd($cmd);
}

##annotate random forest result
my $annotate_matrix=$output_matrix;
$annotate_matrix=~s/\.\w+$/.annotate.txt/;

if($omit_qrf_result_anno eq ""){
	$cmd="paste $matrixqtlresult_hic $output_matrix | cut -f1-5,9,16 > $annotate_matrix\n";
	run_cmd($cmd);
	
}

##delete temprary files
if(($omit_matrixeqtl eq "") and (($matrixqtlresult eq "") or ($matrixqtlresult eq "NULL"))){
		my $tmp="$matrixqtlresult";
		$tmp=~s/cis-qtl/trans-qtl/;
		`unlink $tmp`;
}

if($omit_matrixeqtl_result_anno eq ""){
	`unlink $matrixqtlresult_cor`;
	`unlink $matrixqtlresult_cor_uniq`;
	$matrixqtlresult_cor=~s/\.bed$//;
	`unlink ${matrixqtlresult_cor}.Cpg.bed`;
	`unlink ${matrixqtlresult_cor}.Snp.bed`;
}


if($omit_add_hic eq ""){
	`unlink $matrixqtlresult_hic`;
}


if($omit_add_recomb eq ""){
	`unlink $matrixqtlresult_recomb`;
}


if($omit_cal_recomb_oe eq ""){
	`unlink $matrixqtlresult_recomb_oe`;
}

if($omit_make_input_matrix eq ""){
	`unlink $input_matrix`;
}


if($omit_random_forest eq ""){
	`unlink $output_matrix`;
}

##compare it with MatrixEQTL

##use plink to generate LD independent interval, and find GWAS category enrichment


##########################################All Subroutines###############################################################
sub check_parameter{
	usage() if ( scalar(@ARGV) <= 1 );
	
	if($mode != 1 && $mode != 2){
		print STDERR "Wrong mode number!!\n\n";
		usage();
	}
	
}


sub check_conf_file{
	my $input=shift @_;
	my $dir="";
	my $matrixqtlresult="";
	my $snp_loc="";
	my $gene_loc="";
	my $snp_tab="";
	my $gene_tab="";
	my $cov_tab="";
	my $recomb="";
	my $hic="";
	my $train="";
	
	open(FH,"<$input") or die "can't open configure file $input:$!\n";
	while(<FH>){
		chomp;
		my @f=split "=";
		if($f[0] eq "directory"){
			$dir=$f[1];
		}elsif($f[0] eq "MatrixQtlResult"){
			$matrixqtlresult=$f[1];
		}elsif($f[0] eq "SNP_loc"){
			$snp_loc=$f[1];
		}elsif($f[0] eq "gene_loc"){
			$gene_loc=$f[1];
		}elsif($f[0] eq "SNP_tab"){
			$snp_tab=$f[1];
		}elsif($f[0] eq "gene_tab"){
			$gene_tab=$f[1];
		}elsif($f[0] eq "covar_tab"){
			$cov_tab=$f[1];
		}elsif($f[0] eq "recombination_bw"){
			$recomb=$f[1];
		}elsif($f[0] eq "HiC_KRnorm"){
			$hic=$f[1];
		}elsif($f[0] eq "training_file"){
			$train=$f[1];
		}elsif($f[0] eq "recombination_expect"){
			$recombination_expect=$f[1];
		}
	}
	close(FH);
	
	##check parameter read
	if(($dir eq "") or ($dir eq "NULL")){
		$dir="./";
	}
	if(($recomb eq "") or ($recomb eq "NULL") or ($recombination_expect eq "") or ($recombination_expect eq "NULL")){
		die "recombination rate bigwig file could not be NULL!!\n\n";
	}
	if(($hic eq "") or ($hic eq "NULL")){
		die "HiC file could not be NULL!!\n\n";
	}
	if(($train eq "") or ($train eq "NULL")){
		die "Training file name could not be NULL!!\n\n";
		if(($mode == 1) and ((!-e $train) or (-z $train))){
			die "Training file does not exists or is emptyL!!\n\n";
		}
	}
	if(($matrixqtlresult ne "") and ($matrixqtlresult ne "NULL") and (-e $matrixqtlresult) or (! -z $matrixqtlresult)){
		if(($snp_loc eq "") or ($snp_loc eq "NULL") or ($gene_loc eq "") or ($gene_loc eq "NULL")){
			die "SNP location and Gene/CpG location could not be NULL!!\n\n";
		}
	}else{
		if(($snp_loc eq "") or ($snp_loc eq "NULL") or ($gene_loc eq "") or ($gene_loc eq "NULL") or ($snp_tab eq "") or ($snp_tab eq "NULL") or ($gene_tab eq "") or ($gene_tab eq "NULL")){
			die "SNP location, Gene/CpG location, SNP information and Gene/CpG information could not be NULL!!\n\n";
		}
		if((! -e $matrixqtlresult) or (-z $matrixqtlresult)){
			die "MatrixEQTL result file $matrixqtlresult can not be found or is empty!!\n\n";
		}
	}
	
	return($dir,$matrixqtlresult,$snp_loc,$gene_loc,$snp_tab,$gene_tab,$cov_tab,$recomb, $hic, $train, $recombination_expect);
}



sub run_cmd{
	my $cmd=shift @_;
	print STDERR "$cmd\n";
	system($cmd) == 0 || die "can't execute command $cmd:$!\n";
}

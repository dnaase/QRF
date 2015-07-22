## get unique row of bed files, input bed file should be sorted already.
## author: Yaping Liu  lyping1986@gmail.com

$bed_file=$ARGV[0];

my $use_age = "USAGE: perl uniqBed.pl input_sorted_bed > output_uniq_bed";
if($ARGV[0] eq ""){
	print "$use_age\n";
	exit(1);
}

open(FH,"<$bed_file") or die;
my $chr="";
my $start="";
my $end="";
while(<FH>){
	chomp;
	my $line=$_;
	my @splitin = split "\t",$line;
	if($splitin[0] ne $chr or $splitin[1] ne $start or $splitin[2] ne $end){
		print "$line\n";
		$chr=$splitin[0];
		$start=$splitin[1];
		$end=$splitin[2];	
	}

} 
close(FH);
my $input=$ARGV[0] or die;
my $random=$ARGV[1] or die;  ##the summary file of random permutation
my $observed_data_column=$ARGV[2] || 10;
$observed_data_column--;
my $INFI=100;

my %hash_data=();
open(FH,"<$random") or die;
while(<FH>){
	chomp;
	my @f=split "\t";
	my $key=$f[0];
	$hash_data{$key}=$f[1];
}
close(FH);

open(FH,"<$input") or die;
while(<FH>){
	chomp;
	my @f=split "\t";
	my $len=$f[2]-$f[1];
	my $key="$f[0]:0-$len";
	my $o=$f[$observed_data_column];
	if(exists $hash_data{$key}){
		my $e=$hash_data{$key};
		print "$o\t$e\t";
		if(($o eq "NA") or ($e eq "NA")){
			print "NA\tNA\n";
		}else{
			if($e==0){
				print $INFI."\t";
			}else{
				print $o/$e."\t";
			}
			print $o-$e."\n";
		}
	}
}
close(FH);



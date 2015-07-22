my $qtl=$ARGV[0] or die "no matrix eqtl result file\n";
my $snp_pos=$ARGV[1] or die "no SNP position file\n";
my $cg_pos=$ARGV[2] or die "no CpG/genes position file\n";
my %hash_snp=();
open(F,"<$snp_pos") or die;
while(<F>){
        chomp;
        my @f=split "\t";
        my $l=$f[1]."-".$f[2];
        $hash_snp{$f[0]}=$l;
}
close(F);

my %hash_cg=();
open(F,"<$cg_pos") or die;
while(<F>){
        chomp;
        my @f=split "\t";
        my $l=$f[1]."-".$f[2]."-".$f[3];
        $hash_cg{$f[0]}=$l;
}
close(F);

open(F,"<$qtl") or die "can not open matrix eqtl result file:$!\n";
my $out=$qtl;

$out=~s/\.\w+$/.sig.sameChr.addCor.bed/;

my $out_cpg=$qtl;
$out_cpg=~s/\.\w+$/.sig.sameChr.addCor.Cpg.bed/;
my $out_snp=$qtl;
$out_snp=~s/\.\w+$/.sig.sameChr.addCor.Snp.bed/;

open(OUT,">$out") or die "can not write $out file:$!\n";
open(OUTCG,">$out_cpg") or die "can not write $out_cpg file:$!\n";
open(OUTSNP,">$out_snp") or die "can not write $out_snp file:$!\n";
while(<F>){
        chomp;
        my @f=split "\t";
        if((exists $hash_snp{$f[0]}) and (exists $hash_cg{$f[1]})){
                my @snp_cor=split "-", $hash_snp{$f[0]};
                my @cg_cor=split "-", $hash_cg{$f[1]};
                if($snp_cor[0] eq $cg_cor[0]){
                        my $chr=$snp_cor[0];
                        my $start=$snp_cor[1]<$cg_cor[1]?($snp_cor[1]-1):$cg_cor[1];
                        my $end=$snp_cor[1]>$cg_cor[2]?$snp_cor[1]:$cg_cor[2];
                        print OUT "$chr\t$start\t$end\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\n";
                        print OUTSNP "$chr\t".($snp_cor[1]-1)."\t$snp_cor[1]\t$chr:$start-$end\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\n";
                        print OUTCG "$chr\t$cg_cor[1]\t$cg_cor[2]\t$chr:$start-$end\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\n";
                }
        }

}
close(F);
close(OUT);
close(OUTCG);
close(OUTSNP);
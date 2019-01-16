package align;

use strict;
use warnings;
use File::Basename;
use String::Random;

use FindBin;
use lib $FindBin::Bin;
use Rcall qw ( histogram );

use Exporter;
our @ISA         = qw( Exporter );
our @EXPORT    = qw( &rpms_rpkm_te &BWA_call &to_build &get_unique &sam_sorted_bam &get_hash_alignment &sam_to_bam_bg &sam_count &sam_count_mis &rpms_rpkm &get_fastq_seq &extract_sam );

sub to_build
{
    my ( $toBuildTabP, $log, $newdir ) = @_;

    foreach my    $pairs ( @{ $toBuildTabP } )
    {
        if (    $pairs->[0] == 1 )
        {
            my $sym = $newdir.basename(${$pairs->[1]}).'_symlink.fa';
            symlink( ${$pairs->[1]}, $sym );
            ${$pairs->[1]} = $sym;
            build_index ( ${$pairs->[1]}, $log );
        }
    }
}

sub build_index
{
    my $to_index = shift;
    my $log = shift;
    my $index_log = $to_index.'_index.err';
    `bwa index '$to_index' 2> '$index_log'`;
    print $log "Creating index for $to_index\n";
}

sub get_unique
{
    my ( $sam, $s_uni, $out_prefix, $col_prefix, $details, $report ) = @_;

    my $fout = $col_prefix.'_all_mappers.fastq';
    my $funi = $col_prefix.'_unique_mappers.fastq';
    my $frej = $col_prefix.'_unmapped.fastq';
    
    my $repartition = $out_prefix.'distribution.txt';
    my $png_rep = $out_prefix.'distribution.png';
    my ( %duplicates, %genome_hits) ;

    #alignement to the first reference
    my @return = sam_parse( $sam, $fout, $funi, $frej, $s_uni, \%duplicates, \%genome_hits, $report );
    my $ref_fai = $return[4];
    my $mappers =    $return[5];
    my $mappers_uni = $return[6];
    my $size_mappedHashR = $return[7];

    if ( $details == 1 )
    {
        #print number of duplicates and hits number
        my ($pourcentage, $total) =(0,0);

        $total += $_ foreach values %{$size_mappedHashR};
        open (my $rep, '>'.$repartition) || die "cannot create $repartition $!\n";
        print $rep "size\tnumber\tpercentage\n";
        foreach my $k (sort{$a cmp $b} keys (%{$size_mappedHashR}))
        {
            $pourcentage = 0;
            $pourcentage = $size_mappedHashR->{$k} / $total * 100 unless $total ==0;

            print $rep "$k\t$size_mappedHashR->{$k}\t";
            printf $rep "%.2f\n",$pourcentage;
        }

        histogram($size_mappedHashR, $png_rep, $total);


        my $dup = $out_prefix.'dup_mapnum.txt';
        my $dup_u = $out_prefix .'dup_unique.txt';
        my $dup_r = $out_prefix .'dup_nonmapp.txt';
        open(my $tab,">".$dup) || die "cannot open output txt file\n";
        open(my $tab_r,">".$dup_r) || die "cannot open output txt file\n";
        open(my $tab_u,">".$dup_u) || die "cannot open output txt file\n";
        print $tab "sequence\tcount\tmapnum\n";
        print $tab_u "sequence\tcount\n";
        print $tab_r "sequence\tcount\n";
        foreach my $k (sort {$duplicates{$b} <=> $duplicates{$a}}keys %duplicates)
        {
            $duplicates{$k} = 0 unless exists($duplicates{$k});
            $genome_hits{$k} = 0 unless exists($genome_hits{$k});
            if ($genome_hits{$k} != 0) { print $tab $k."\t".$duplicates{$k}."\t".$genome_hits{$k}."\n"; }
            else {print $tab_r $k."\t".$duplicates{$k}."\n";}
            if ($genome_hits{$k} == 1) { print $tab_u $k."\t".$duplicates{$k}."\n"; }
        }
    close $dup; close $dup_r; close $dup_u;
    }
    return ( $ref_fai, $mappers, $mappers_uni );
}

sub sam_parse
{
    my ( $sam, $fastq_accepted, $fastq_accepted_unique, $fastq_rejected, $sam_unique, $duplicate_hashR, $best_hit_number_hashR, $report ) = @_ ;
    my ($reads, $mappers, $mappersUnique, @garbage, %size_num, %size_num_spe, %number, %numberSens, %numberReverse, %unique_number, %numberNM, %numberM, %size);
    $mappers = $mappersUnique = $reads = 0;

    open my $fic, '<', $sam || die "cannot open $sam $!\n";
    open my $accepted, '>', $fastq_accepted || die "cannot create $fastq_accepted $! \n";
    open my $unique, '>', $fastq_accepted_unique || die "cannot create $fastq_accepted_unique $! \n";
    open my $rejected, '>', $fastq_rejected || die "cannot create $fastq_rejected $! \n";
    open my $sam_uni, '>', $sam_unique || die "cannot create $sam_unique $! \n";
    my $sequence = '';
    while(<$fic>)
    {
        chomp $_;
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                $size{$1} = $2;
                $unique_number{$1} = 0;
                $number{$1} = 0;
                $numberNM{$1} = 0;
                $numberM{$1} = 0;
            }
            print $sam_uni $_."\n";
            next;
        }
        $reads++;
        my @line = split (/\t/,$_);
        $sequence = $line[9];
        if ($line[1] & 16)
        {
            $sequence =reverse($sequence);
            $sequence =~ tr/atgcuATGCU/tacgaTACGA/;
        }
        if ($line[1] == 16 || $line[1] == 0)
        {
            my $len = length($sequence);
            $size_num{$len} ++;
            $size_num_spe{$line[2]}{$len}++;
            $mappers ++;

            ${$best_hit_number_hashR}{$sequence} = $1    if    ($line[13] =~ /X0:i:(\d*)/ ||    $line[14] =~/X0:i:(\d*)/ );
            ${$duplicate_hashR}{$sequence}++;
            $number{$line[2]}++;
            $numberSens{$line[2]}++ if $line[1] == 0 ;
            $numberReverse{$line[2]}++ if $line[1] == 16 ;
            print $accepted "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";

            if ($line[11] eq "XT:A:U")
            {
                $unique_number{$line[2]}++;
                $mappersUnique++;
                print $unique "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
                print $sam_uni $_."\n";
            }
            if ($_ =~ /.*XM:i:(\d+).*/)
            {
                if ($1 == 0){$numberNM{$line[2]}++;}else{$numberM{$line[2]}++;}
            }
        }
        else
        {
            ${$best_hit_number_hashR}{$sequence} = 0;
            ${$duplicate_hashR}{$sequence}++;
            print $rejected "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
        }
    }
    close $fic; close $accepted; close $unique; close $rejected; close $sam_uni;

    print $report "Parsing $sam file\n";
    print $report "\treads: $reads\n";
    print $report "\tmappers: $mappers\n";
    print $report "\tunique mappers: $mappersUnique\n";
    print $report "-----------------------------\n";
    return (\%number, \%unique_number, \%numberSens, \%numberReverse, \%size, $mappers, $mappersUnique, \%size_num, \%size_num_spe, \%numberNM, \%numberM );
}

sub get_hash_alignment
{
    my ($index, $mismatches, $accept, $reject, $outA, $outR, $fastq, $number_of_cpus, $name, $sam, $report, $fai_f) = @_ ;
    my ($reads, $mappers, $unmapped) = (0,0,0);
    my $accep_unique;
    BWA_call ( $index, $fastq, $sam, $mismatches, $number_of_cpus, $report );
    
    open my $fic, '<', $sam || die "cannot open $sam $!\n";
    open my $accepted, '>', $outA || die "cannot open $outA\n"    if $accept == 1;
    open my $rejected, '>', $outR || die "cannot open $outR\n"    if $reject == 1;
    open my $fai, '>', $fai_f || die "cannot open $fai_f\n" if $fai_f;
    #if ($name eq "snRNAs") {
    #    open ( $accep_unique, ">".$1."-unique.fastq") if $outR =~ /(.*)\.fastq/;
    #}
    my $sequence = '';
    while(<$fic>)
    {
        chomp $_;
        if( $_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($fai_f && $_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                print $fai $1."\t".$2."\n";
            }
            next;
        }
        $reads++;
        my @line = split (/\t/,$_);
        $sequence = $line[9];
        if ($line[1] & 16)
        {
            $sequence =reverse($sequence);
            $sequence =~ tr/atgcuATGCU/tacgaTACGA/;
        }
        if ($line[1] & 16 || $line[1] == 0)
        {
            $mappers ++;
            if ($accept == 1 )
            {
                print $accepted "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n";
     #         print $accep_unique "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n" if ($name eq "snRNAs" && $line[11] eq "XT:A:U");
            }
        }
        else
        {
            print $rejected "\@".$line[0]."\n".$sequence."\n+\n".$line[10]."\n" if $reject == 1;
            $unmapped++;
        }
    }
 # close $accep_unique if ($name eq "bonafide_reads");
    close $fic;
    close $accepted    if $accept == 1;
    close $rejected if $reject ==1;
    close $fai if $fai_f;
    print $report "\treads: $reads\n";
    print $report "\tmappers: $mappers\n";
    print $report "\tunmapped: $unmapped\n";
    print $report "-----------------------------\n";
    return ($mappers, $unmapped);
}

sub sam_count
{
    my $sam = shift;
    my ( %number, %size );
    
    open    my $fic, '<', $sam || die "cannot open $sam file $!\n";
    while(<$fic>)
    {
        chomp $_;
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                $size{$1} = $2;
                $number{$1} = 0;
            }
        }
        else
        {
            my @line = split (/\t/,$_);
             if ( $line[1]    & 16 || $line[1] == 0 )
             {
                 $number{$line[2]}++;
             }
          }
    }
    close $fic;
    return ( \%number, \%size );
}

sub sam_count_mis
{

    my $sam = shift;
    my ( %number, %numberNM, %numberM, %size);
    
    open    my $fic, '<', $sam || die "cannot open $sam file $!\n";
    while(<$fic>)
    {
        chomp $_;
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            if ($_ =~ /\@SQ\tSN:(.*)\tLN:(\d*)/)
            {
                $size{$1} = $2;
                $number{$1} = [0,0,0,0,0,0,0];
                $numberNM{$1} = [0,0,0,0,0,0,0];
                $numberM{$1} = [0,0,0,0,0,0,0];
            }
        }
        else
        {
            my @line = split (/\t/,$_);
            my @seq = split //, $line[9];
            if ( $line[1]    == 16 || $line[1] == 0 )
            {
                $number{ $line[2] }->[0]++;
                if ($line[1]    == 0)
                {
                    $number{$line[2]}->[1]++;
                    $number{$line[2]}->[3]++ if $seq[0] eq 'T';
                    $number{$line[2]}->[5]++ if $seq[9] eq 'A';
                }
                else
                {
                    $number{$line[2]}->[2]++;
                    $number{$line[2]}->[4]++ if $seq[9] eq 'A';
                    $number{$line[2]}->[6]++ if $seq[0] eq 'T';
                }
                 if ($_ =~ /.*XM:i:(\d+).*/)
                 {
                    if ( $1 == 0 )
                    {
                        $numberNM{$line[2]}->[0]++;
                        if ($line[1]    == 0)
                        {
                            $numberNM{$line[2]}->[1]++;
                            $numberNM{$line[2]}->[3]++ if $seq[0] eq 'T';
                            $numberNM{$line[2]}->[5]++ if $seq[9] eq 'A';
                        }
                        else
                        {
                            $numberNM{$line[2]}->[2]++;
                            $numberNM{$line[2]}->[4]++ if $seq[9] eq 'A';
                            $numberNM{$line[2]}->[6]++ if $seq[0] eq 'T';
                        }
                    }
                    else
                    {
                        $numberM{$line[2]}->[0]++;
                        if ($line[1]    == 0)
                        {
                            $numberM{$line[2]}->[1]++;
                            $numberM{$line[2]}->[3]++ if $seq[0] eq 'T';
                            $numberM{$line[2]}->[5]++ if $seq[9] eq 'A';
                        }
                        else
                        {
                            $numberM{$line[2]}->[2]++;
                            $numberM{$line[2]}->[4]++ if $seq[9] eq 'A';
                            $numberM{$line[2]}->[6]++ if $seq[0] eq 'T';
                        }
                    }
                }
            }
        }
    }
    return (\%number, \%size, \%numberNM, \%numberM );
}

sub rpms_rpkm_te
{
    my ( $counthashR, $sizehashR, $mapped, $out_file, $piRNA_number, $miRNA_number, $bonafide_number ) =@_;
    open(my $out, ">".$out_file) || die "cannot open normalized file $! \n";
    print $out "ID\treads counts\tRPKM";
    print $out "\tper million of piRNAs" if ($piRNA_number != 0);
    print $out "\tper million of miRNAs" if ($miRNA_number != 0);
    print $out "\tper million of bonafide reads" if ($bonafide_number != 0);
    print $out "\tsense reads counts\treverse reads counts";
    print $out "\t% of sense 1U\t% of sense 10A\t% of reverse 1U\t% of reverse 10A\n";
    foreach my $k    ( sort keys %{$counthashR} )
    {
        my ($rpkm, $pirna, $mirna, $bonafide) = (0,0,0,0);
        
        $rpkm = ( $counthashR->{$k}->[0] * 1000000000) / ( $sizehashR->{$k} * $mapped) if ( $sizehashR->{$k} * $mapped) != 0 ;
        print $out $k."\t".$counthashR->{$k}->[0]."\t"; printf $out "%.2f",$rpkm;

        if ($piRNA_number != 0 )
        {
            $pirna = ( $counthashR->{$k}->[0]    * 1000000) / $piRNA_number;
            printf $out "\t%.2f",$pirna;
        }
        if ($miRNA_number != 0 )
        {
            $mirna = ( $counthashR->{$k}->[0]    * 1000000) / $miRNA_number;
            printf $out "\t%.2f",$mirna;
        }
        if ($bonafide_number != 0 )
        {
            $bonafide = ( $counthashR->{$k}->[0]    * 1000000) / $bonafide_number;
            printf $out "\t%.2f",$bonafide;
        }

        print $out "\t".$counthashR->{$k}->[1]."\t".$counthashR->{$k}->[2] ;
        my $S1U = 0;
        $S1U = $counthashR->{$k}->[3] / $counthashR->{$k}->[1] * 100 if $counthashR->{$k}->[1] != 0;
        my $R1U = 0;
        $R1U = $counthashR->{$k}->[6] / $counthashR->{$k}->[2] * 100 if $counthashR->{$k}->[2] != 0;
        my $S10A = 0;
        $S10A = $counthashR->{$k}->[5] / $counthashR->{$k}->[1] * 100 if $counthashR->{$k}->[1] != 0;
        my $R10A = 0;
        $R10A = $counthashR->{$k}->[4] / $counthashR->{$k}->[2] * 100 if $counthashR->{$k}->[2] != 0;
        print $out "\t".$S1U."\t".$S10A."\t".$R1U."\t".$R10A;

        print $out "\n";
    }
    close $out;
}


sub sam_to_bam_bg
{
    my ( $sam, $scale, $number_of_cpus ) = @_;
    my ( $bam_sorted, $bedgraphM, $bedgraphP, $view_err, $sort_err ) = ( '', '', '', '', '' );
    if ( $sam =~ /(.*?).sam$/ )
    {
        $bam_sorted = $1.'_sorted.bam';
        $bedgraphP= $1.'_plus.bedgraph';
        $bedgraphM = $1.'_minus.bedgraph';
        $view_err = $1.'_view.err';
        $sort_err = $1.'_sort.err';
    }
    `samtools view -Shb     --threads $number_of_cpus '$sam' 2> '$view_err' | samtools sort    -O BAM --threads $number_of_cpus    /dev/stdin 2> '$sort_err'    > '$bam_sorted'`;
    `bedtools genomecov -scale $scale -strand + -bga -ibam '$bam_sorted' > '$bedgraphP'`;
    `bedtools genomecov -scale $scale -strand - -bga -ibam '$bam_sorted' > '$bedgraphM'`;
}

sub sam_sorted_bam
{
    my ( $sam, $number_of_cpus ) = @_;
    my ( $bam_sorted, $view_err, $sort_err ) = ( '', '', '' );
    if ( $sam =~ /(.*?).sam$/ )
    {
        $bam_sorted = $1.'_sorted.bam';
        $view_err = $1.'_view.err';
        $sort_err = $1.'_sort.err';

    }
    `samtools view -Shb     --threads $number_of_cpus '$sam' 2> '$view_err' | samtools sort    -O BAM --threads $number_of_cpus    /dev/stdin    2> '$sort_err'    > '$bam_sorted'`;
}

sub BWA_call
{
    my ( $index, $fastq, $sam, $mismatches, $number_of_cpus, $report ) = @_;
    my ( $aln_err, $samse_err, $seq_num ) = ( $sam.'_aln.err', $sam.'_samse.err', 0 );
    print $report "-----------------------------\n";
    print $report "bwa aln -t $number_of_cpus -n $mismatches '$index' '$fastq' 2> '$aln_err' | bwa samse $index /dev/stdin '$fastq' 2> '$samse_err' > '$sam'\n";
    `bwa aln -t $number_of_cpus -n $mismatches '$index' '$fastq' 2> '$aln_err' | bwa samse $index /dev/stdin '$fastq' 2> '$samse_err' > '$sam' `;
}

sub rpms_rpkm
{
    my ( $counthashR, $sizehashR, $mapped, $out_file, $piRNA_number, $miRNA_number, $bonafide_number ) =@_;
    open(my $out, ">".$out_file) || die "cannot open normalized file $! \n";
    print $out "ID\treads counts\tRPKM";
    print $out "\tper million of piRNAs" if ($piRNA_number != 0);
    print $out "\tper million of miRNAs" if ($miRNA_number != 0);
    print $out "\tper million of bonafide reads" if ($bonafide_number != 0);
    print $out "\n";
    foreach my $k    ( sort keys %{$counthashR} )
    {
        my ($rpkm, $pirna, $mirna, $bonafide) = (0,0,0,0);
        
        $rpkm = ( $counthashR->{$k} * 1000000000) / ( $sizehashR->{$k} * $mapped) if ( $sizehashR->{$k} * $mapped) != 0 ;
        print $out $k."\t".$counthashR->{$k}."\t"; printf $out "%.2f",$rpkm;
        
        if ($piRNA_number != 0 )
        {
            $pirna = ( $counthashR->{$k}    * 1000000) / $piRNA_number;
            printf $out "\t%.2f",$pirna;
        }
        if ($miRNA_number != 0 )
        {
            $mirna = ( $counthashR->{$k}    * 1000000) / $miRNA_number;
            printf $out "\t%.2f",$mirna;
        }
        if ($bonafide_number != 0 )
        {
            $bonafide = ( $counthashR->{$k}    * 1000000) / $bonafide_number;
            printf $out "\t%.2f",$bonafide;
        }
        print $out "\n";
     }
    close $out;
}

sub extract_sam
{
    my ( $hashRef, $sam_in, $sam_out, $sam_uni_out, $fastq_out, $fastq_uni_out ) = @_;
    
    open my $s_in, '<', $sam_in || die "cannot open $sam_in file $!\n";

    open my $f_out, '>', $fastq_out || die "cannot create $fastq_out $!\n";
    open my $f_uni_out, '>', $fastq_uni_out || die "cannot create $fastq_uni_out $!\n";
    
    open my $s_out, '>', $sam_out || die "cannot create $sam_out file $!\n" if defined ($hashRef);
    open my $s_uni_out, '>', $sam_uni_out || die "cannot create $sam_uni_out file $!\n";

    my $sequence = '';
    while(<$s_in>)
    {
        if ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
        {
            print $s_out $_ if defined ($hashRef);
            print $s_uni_out $_;
            next;
        }
        my @line = split (/\t/,$_);
        $sequence = $line[0];
        if ( (! defined ($hashRef) )|| (    exists $hashRef->{$sequence}    &&    $hashRef->{$sequence} == 1 ) )
        {
            my $arn    =    $line[9];
            if ($line[1] & 16)
            {
                $arn =reverse($arn);
                $arn =~ tr/atgcuATGCU/tacgaTACGA/;
            }

            if    ( ( $line[1] == 16 || $line[1] == 0 ) )
            {
                print $f_out "\@".$line[0]."\n".$arn."\n+\n".$line[10]."\n" ;
                print $s_out $_ if defined ($hashRef);
                if ( $line[11] eq "XT:A:U" )
                {
                    print $f_uni_out "\@".$line[0]."\n".$arn."\n+\n".$line[10]."\n" ;
                    print $s_uni_out $_ ;
                }
            }
        }
    }
    close $s_in; close $s_out if defined ($hashRef);
    close $s_uni_out; close $f_out; close $f_uni_out;
}

sub get_fastq_seq
{
    my $fastq = shift;
    my %hash; my $cmp = 0;

    open my $fic, '<', $fastq || die "cannot open input file $! \n";
    while(<$fic>)
    {
        chomp $_;
        $cmp++;
        if ($cmp % 4 == 1)
        {
            die "file do not contain a @ at line $cmp\n" unless ($_ =~ /^\@/ );
            if ($_ =~ /^\@(.*)\s.*/) { $hash{$1} = 1;}
            elsif ($_ =~ /^\@(.*)/) { $hash{$1} = 1;}
        }
        elsif ($cmp % 4 == 3 )
        {
            die "file do not contain a + at line $cmp\n" unless $_ =~ /^\+/;
        }
    }
    close $fic;
    return \%hash;
}

1;

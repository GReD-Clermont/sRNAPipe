package sRNAPipe::resize;

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use sRNAPipe::Rcall qw ( histogram );

use Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( &size_distribution );

sub size_distribution
{
    my ( $fastq, $fastq_out, $dir, $min, $max ) = @_;

    my ( %fragments_size, %duplicates ) ;
    my $num = size($min, $max, $fastq, $fastq_out, \%fragments_size, \%duplicates);

    my $png = $dir.'histogram.png';
    histogram(\%fragments_size, $png, $num);
    
    my $size = $dir.'reads_size.txt';
    

    my $pourcentage;
    open    my $o, '>', $size || die "cannot open $size $!\n";
    print $o "size\tnumber\tpercentage\n";
    foreach my    $k (sort { $a <=> $b } keys %fragments_size )
    {
        $pourcentage = $fragments_size{$k} / $num * 100;
        
        print $o "$k\t$fragments_size{$k}\t";
        printf $o "%.2f\n",$pourcentage;
    }
    close $o;
    
    my $dup = $dir.'duplicates.txt' ;
    open $o, '>', $dup || die "cannot open $size $!\n";
    print $o "size\tnumber\n";
    foreach my    $k (sort { $duplicates{$b} <=> $duplicates{$a} } keys %duplicates )
    {
        print $o "$k\t$duplicates{$k}\n";
    }
    close $o;
}

sub size
{
    my ($min, $max, $in_file, $out_file, $sizeHashR, $duplicateHashR) = @_;
    my ($numreads, $size, $cmp, $ok, $line) = (0, 0, 0, 0);
    my @fastq;
    open (my $in, $in_file) || die "cannot open $in_file $!\n";
    open (my $out, ">".$out_file)    || die "cannot create $out_file $!\n";
    while(<$in>)
    {
        chomp $_;
        $cmp++; $line++;
        if ($cmp == 1)
        {
            die "file do not contain a @ at line $line\n" unless ($_ =~ /^\@/ );
            $ok = 0; @fastq = ();
            push(@fastq,$_);
        }
        elsif ($cmp == 2)
        {
            #die "unrecognized symbol at line $line\n" unless ($_ =~ /[atcgATCGnN]+/ || $_ =~ /^$/ );
            push(@fastq,$_);
            $size = length($_);
            if ($size >= $min && $size <= $max)
            {
                $numreads++;
                ${$sizeHashR}{$size}+=1;
                ${$duplicateHashR}{$_}+=1 if (defined($duplicateHashR));
                $ok = 1;
            }
        }
        elsif ($cmp == 3 )
        {
            die "file do not contain a + at line $line\n" unless $_ =~ /^\+/;
            push(@fastq,$_);
        }
        elsif ($cmp == 4 )
        {
            push(@fastq,$_);
            $cmp = 0;
            if ($ok == 1)
            {
                foreach my $t (@fastq)
                {
                    print $out $t."\n";
                }
            }
        }
    }
    close $in; close $out;
    return $numreads;
}

1;

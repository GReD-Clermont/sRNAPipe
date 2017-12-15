package ppp;

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use Rcall qw ( histogram );
use Math::CDF;

use Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( &ping_pong_partners );

sub ping_pong_partners
{
  my ( $TE_fai, $sam, $dir, $max ) = @_;
  my ( $hashRef, $dupRef, $hasPpp ) = count_mapped ( $TE_fai, $sam );
  my ( %num_per_overlap_size, $overlap_number, $reverseR, $begRev, $endRev, $sensR, $begSens, $endSens, $snum, $rnum, $overlap );
  my ( $SP, $AP, $SN, $AN, $txt );
  my $flag = 0;
  my @distri_overlap = (); my @overlaps_names = ();

  open my $ppp_f, '>', $dir."ppp.txt" || die "cannot create ppp.txt $!\n";
  foreach my $k ( sort keys %{$hashRef} )
  {
    my $v = $hashRef->{$k};
    my $TE_dir = $dir.$k.'/';

    %num_per_overlap_size = (); $overlap_number = 0;
    $flag = 0;
    for ( my $i = 0; $i <= $#{$v->[1]} ; $i++ )
    {
      $reverseR = ${$v->[1]}[$i] ;
      $begRev = $reverseR->[0];
      $endRev = $begRev + length($reverseR->[1]) - 1;

      my $revR = reverse($reverseR->[1]);
      $revR =~ tr/atgcuATGCU/tacgaTACGA/;

      for ( my $j = 0; $j <= $#{$v->[0]}; $j++ )
      {
        $sensR = ${$v->[0]}[$j];
        $begSens = $sensR->[0];
        $endSens = $begSens + length($sensR->[1]) - 1;

        if ( $begSens <= $endRev && $endSens > $endRev )
        {
          $flag = 1;
          mkdir $TE_dir;
          open  $txt, '>', $TE_dir.'overlap_size.txt' || die "cannot open repartition\n";

          $overlap = $endRev - $begSens + 1;
          $snum =  $dupRef->{$sensR->[0].$sensR->[1].$sensR->[2].$sensR->[3]};
          $rnum = $dupRef->{$reverseR->[0].$reverseR->[1].$reverseR->[2].$reverseR->[3]};

          if ( $overlap == 10 )
          {
            $hasPpp->{ $sensR->[0].$sensR->[1].$sensR->[2].$sensR->[3] } = 1;
            $hasPpp->{ $reverseR->[0].$reverseR->[1].$reverseR->[2].$reverseR->[3] } = 1;
          }
          next if $overlap > $max;
          if ( $snum < $rnum )
          {
            $num_per_overlap_size{$overlap} += $snum;
            $overlap_number += $snum;
          }
          else
          {
            $num_per_overlap_size{$overlap} += $rnum ;
            $overlap_number += $rnum ;
          }
        }
      }
    }
    if ( $max != 0 )
    {
      my @overlaps = ();
      push @overlaps_names, $k;
      for my $i (1..$max)
      {
        $num_per_overlap_size{$i} = 0 unless exists( $num_per_overlap_size{$i} );
        push @overlaps, $num_per_overlap_size{$i};
      }
      push @distri_overlap, \@overlaps;
    }

    if ( $flag == 1 )
    {
      open  $AP, '>', $TE_dir."antisensPPP.txt" || die "cannot create antisensPPP\n";
      open  $AN, '>', $TE_dir."antisens.txt"  || die "cannot create antisens\n";
      for ( my $i = 0; $i <= $#{$v->[1]} ; $i++ )
      {
        $reverseR = ${$v->[1]}[$i] ;
        my $revR = reverse($reverseR->[1]);
        $revR =~ tr/atgcuATGCU/tacgaTACGA/;
        $rnum = $dupRef->{$reverseR->[0].$reverseR->[1].$reverseR->[2].$reverseR->[3]};
        if ( $hasPpp->{ $reverseR->[0].$reverseR->[1].$reverseR->[2].$reverseR->[3] } == 1 )
        {
          print $AP ">$reverseR->[0]|$reverseR->[2]|$reverseR->[3]|$rnum\n$revR\n";
        }
        else
        {
          print $AN ">$reverseR->[0]|$reverseR->[2]|$reverseR->[3]|$rnum\n$revR\n";
        }
      }
      close $AP; close $AN;

      open  $SP, '>', $TE_dir."sensPPP.txt" || die "cannot create sensPPP\n";
      open  $SN, '>', $TE_dir."sens.txt"  || die "cannot create sens\n";
      for ( my $j = 0; $j <= $#{$v->[0]}; $j++ )
      {
        $sensR = ${$v->[0]}[$j];
        $snum =  $dupRef->{$sensR->[0].$sensR->[1].$sensR->[2].$sensR->[3]};
        if ( $hasPpp->{ $sensR->[0].$sensR->[1].$sensR->[2].$sensR->[3] } == 1 )
        {
          print $SP ">$sensR->[0]|$sensR->[2]|$sensR->[3]|$snum\n$sensR->[1]\n";
        }
        else
        {
          print $SN ">$sensR->[0]|$sensR->[2]|$sensR->[3]|$snum\n$sensR->[1]\n";
        }
      }
      close $SP; close $SN;

      my $histo_png = $TE_dir.'histogram.png';
      histogram( \%num_per_overlap_size, $histo_png, $overlap_number );
      print $txt "size\tnumber\tpercentage of the total overlap number\n";
      foreach my $k ( sort {$a <=> $b} keys %num_per_overlap_size )
      {
        my $percentage = 0;
        $percentage = $num_per_overlap_size{$k} * 100 / $overlap_number unless $overlap_number == 0;
        print $txt "$k\t$num_per_overlap_size{$k}\t"; printf $txt "%.2f\n",$percentage;
      }
      close $txt;
    }
  }

  foreach my $tabP (  @distri_overlap )
  {
    my $sum = sum($tabP);
    my $ten = $tabP->[9];
    my $mean = mean($tabP);
    my $std = standard_deviation($tabP, $mean);
    my $zsc = z_significance($ten, $mean, $std);
    my $name = shift @overlaps_names;
    my $prob = 'NA';
    $prob =  1 - &Math::CDF::pnorm( $zsc ) if $zsc ne 'NA';
    print $ppp_f (join ("\t", $name, $sum, $ten, $mean, $std, $zsc, $prob ),"\n" );
  }
  close $ppp_f;
}

sub count_mapped
{
  my ( $fai, $in_file ) = @_;
  my ( %mapped, %dup, %has_ppp );

  open my $f, '<', $fai || die "cannot open $fai $! \n";
  while(<$f>)
  {
    if ($_ =~ /(.*)\t(\d+)\n/)
    {
      $mapped{$1} = [];
      $mapped{$1}->[0] = []; $mapped{$1}->[1] = [];
    }
  }
  close $f;

  open my $infile, "samtools view  $in_file |"|| die "cannot open input file $! \n";
  while(<$infile>)
  {
    unless ($_ =~ /^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$/ || $_ =~ /^\@CO\t.*/ )
    {
      my @line = split (/\t/,$_);
      if ($line[1] == 0)
      {
        unless ( exists ($dup{$line[3].$line[9].$line[1].$line[2]}) )
        {
          push @{$mapped{$line[2]}->[0]} , [$line[3], $line[9], $line[1],  $line[2]];
          $has_ppp {$line[3].$line[9].$line[1].$line[2]} = 0;
        }
        $dup{$line[3].$line[9].$line[1].$line[2]}+=1;
      }
      elsif ($line[1] == 16)
      {
        unless ( exists ($dup{$line[3].$line[9].$line[1].$line[2]}) )
        {
          push @{$mapped{$line[2]}->[1]} , [$line[3], $line[9], $line[1],  $line[2]];
          $has_ppp{$line[3].$line[9].$line[1].$line[2]} = 0;
        }
        $dup{$line[3].$line[9].$line[1].$line[2]}+=1
      }
    }
  }
  close $infile;
  return (\%mapped, \%dup, \%has_ppp );
}

sub sum
{
  my $arrayref = shift;
  my $result = 0;
  foreach (@$arrayref) {$result += $_}
  return $result;
}

sub mean
{
  my $arrayref = shift;
  my $result;
  foreach (@$arrayref) {$result += $_}
  return $result / scalar(@$arrayref);
}

sub standard_deviation
{
  my ($arrayref, $mean) =  @_;
  return sqrt ( mean ( [map $_**2 , @$arrayref ]) - ($mean**2));
}

sub z_significance
{
  my ($ten, $mean, $std) = @_;
  my $z = 'NA';
  $z = (($ten - $mean) / $std) if $std != 0;
  return $z;
}

1;


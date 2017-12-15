package subgroups;

use strict;
use warnings;
use Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( &subgroups );

use POSIX;
use FindBin;
use lib $FindBin::Bin;
use align qw ( get_hash_alignment );

sub subgroups
{
	my ($fin, $dir, $mis, $mis_TE, $proc, $tRNAs, $rRNAs, $snRNAs, $miRNAs, $transcripts, $TE, $min_si, $max_si, $min_pi, $max_pi, $report ) = @_;
	my (@files, $sum, $pie, $repar, %ismapped, %isjunk, %repartition, @junk_ref, @all_ref );

	srand();
	print $report "----------------------------\n";
	print $report "Create subgroups:\nfastq_in: $fin\ndirectory_out: $dir\nmismatches: $mis\nmismatches TE: $mis_TE\n";

	mkdir $dir;
	$dir = $dir.'/' unless $dir =~ /(.*)\/$/;

	my $accept_miRNas = $dir.'miRNAs.fastq'; 
	my $reject_miRNAs = $dir.'miRNAs_rejected.fastq';
	my $sam_miRNAs = $dir.'miRNAs.sam'; 
	my @tmp = get_hash_alignment($miRNAs, $mis, 1, 1, $accept_miRNas, $reject_miRNAs, $fin, $proc, 'miRNAs',$sam_miRNAs, $report);
	my $mi = $tmp[0];
	$repartition{'miRNAs'} = $mi;

	my $sam = new String::Random;
	$sam = $sam->randpattern("CCcccccc");
	my $reject_rRNAs = $dir.'rRNAs_rejected.fastq';
	@tmp = get_hash_alignment($rRNAs, $mis, 0, 1, 'NA', $reject_rRNAs, $reject_miRNAs, $proc, 'rRNAs', $sam, $report);
	$repartition{'rRNAs'} = $tmp[0];
	unlink $sam, $sam.'_aln.err', $sam.'_samse.err';

	$sam = new String::Random;
	$sam = $sam->randpattern("CCcccccc");
	my $reject_tRNAs = $dir.'tRNAs_rejected.fastq';
	@tmp = get_hash_alignment($tRNAs, $mis, 0, 1, 'NA', $reject_tRNAs, $reject_rRNAs, $proc, 'tRNAs', $sam, $report);
	$repartition{'tRNAs'} = $tmp[0];
	unlink $sam, $sam.'_aln.err', $sam.'_samse.err';

	$sam = new String::Random;
	$sam = $sam->randpattern("CCcccccc");
	my $bonafide = $dir.'bonafide_reads.fastq';
	@tmp = get_hash_alignment($snRNAs, $mis, 0, 1, 'NA', $bonafide, $reject_tRNAs, $proc, 'snRNAs', $sam, $report);
	$repartition{'snRNAs'} = $tmp[0];
	my $bo = $tmp[1];
	unlink $sam, $sam.'_aln.err', $sam.'_samse.err';

	my $sam_transcripts = $dir.'transcripts.sam'; 
	my $reject_transcripts = $dir.'rejected_transcripts.fastq';
	@tmp = get_hash_alignment($transcripts, $mis, 0, 1, 'NA', $reject_transcripts, $bonafide, $proc, 'transcripts', $sam_transcripts, $report, $dir.'transcripts.fai');
	$repartition{'transcripts'} = $tmp[0];

	
	my $sam_TEs = $dir.'TEs.sam';
	my $reject_TEs = $dir.'rejected.fastq';
	@tmp = get_hash_alignment($TE, $mis_TE, 0, 1, 'NA', $reject_TEs, $reject_transcripts, $proc, 'TEs', $sam_TEs, $report, $dir.'TEs.fai' );
	$repartition{'TEs'} = $tmp[0] ; $repartition{'others'} = $tmp[1];
	unlink $sam, $sam.'_aln.err', $sam.'_samse.err';
	unlink $reject_transcripts;
	unlink $reject_rRNAs;
	unlink $reject_miRNAs;
	unlink $reject_tRNAs;

	#create repartition
	my $pi = fastqSubgroups($bonafide, $dir, $min_si, $max_si, $min_pi, $max_pi );

	open (my $re, '>'.$dir.'repartition.txt') || die "cannot open $dir repartition.txt $!\n";
	print $re "type\tnumber\tpercentage\n";
	$sum += $_ foreach values %repartition;
	foreach my $k  ( sort keys  %repartition ) 
	{
		my $prct = 0;
		$prct = $repartition{$k} / $sum * 100 if $sum != 0;
		print $re "$k\t$repartition{$k}\t"; printf $re "%.2f\n",$prct;
	}
	return ( $bo, $mi, $pi);
}

sub fastqSubgroups
{
  my ( $fastq, $output_directory, $min_si, $max_si, $min_pi, $max_pi ) = @_;
  my $fastq_siRNA = $output_directory."siRNAs.fastq";
  my $fastq_piRNA = $output_directory."piRNAs.fastq";

  open my $fic, '<', $fastq || die "cannot open input file $! \n";
  open my $si, '>', $fastq_siRNA || die "cannot open siRNA.fastq $! \n";
  open my $pi, '>', $fastq_piRNA || die "cannot open piRNA.fastq $! \n";
  
  my ($length, $cmp, $type, $siRNA_number, $miRNA_h_number, $piRNA_number, $not_pi_number) = (0,0,0,0,0,0,0);
  my (@fastq) =(); my $seq_name;
  my $out;
  while(<$fic>)
  {
    chomp $_;
    $cmp++;
    if ($cmp == 1)
    {
      die "file do not contain a @ at line $cmp\n" unless ($_ =~ /^\@/ );
      $type = 0; @fastq = ();
      if ($_ =~ /^\@(.*)\s.*/) { $seq_name = $1;}
      elsif ($_ =~ /^\@(.*)/) {$seq_name = $1;}
      push(@fastq,$_);
    }
    elsif ($cmp == 2)
    {
      #die "unrecognized symbol at line $cmp\n" unless $_ =~ /[atcgATCGnN]+/;
      push(@fastq,$_);
      $length = length($_);
      $type = 0;
      if ( $length >= $min_si  && $length <= $max_si )
      {
         $type = 2;
         $siRNA_number++;
      }
      if ($length >= $min_pi  && $length <= $max_pi )
      {
        $type += 4;
        $piRNA_number++;
      }      
    }
    elsif ($cmp == 3 )
    {
      die "file do not contain a + at line $cmp\n" unless $_ =~ /^\+/;
      push(@fastq,$_);
    }
    elsif ($cmp == 4 )
    {
      push(@fastq,$_);
      $cmp = 0;
      if ($type != 0)
      {
        if ($type & 4 ) { foreach my $t (@fastq) { print $pi $t."\n";} }
        if ($type & 2 ) { foreach my $t (@fastq) { print $si $t."\n";} }        
      }
    }
  }
  close $fic;
  close $si; close $pi;
  return ($piRNA_number);
}

1;

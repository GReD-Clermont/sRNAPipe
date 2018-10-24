#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use File::Basename;
use File::Copy::Recursive qw( dircopy );
use POSIX;
use FindBin;
use lib $FindBin::Bin;
use resize qw ( size_distribution );
use subgroups qw (subgroups );
use ppp qw ( ping_pong_partners );
use Rcall qw (pie_chart bg_to_png );
use align qw ( to_build get_unique sam_count sam_count_mis sam_sorted_bam rpms_rpkm rpms_rpkm_te BWA_call get_fastq_seq extract_sam sam_to_bam_bg );
use html qw ( main_page details_pages menu_page ppp_page );
use File::Copy;

my ( @fastq, @fastq_n, $dir, $min, $max, $mis, $misTE, $help, $Pcheck, $mapnumf, $html_out);
my ( $ref, $tRNAs, $rRNAs, $snRNAs, $miRNAs, $transcripts, $TE );
my ( $si_min, $si_max, $pi_min, $pi_max );
my ( $build_index, $build_tRNAs, $build_rRNAs, $build_snRNAs, $build_miRNAs, $build_transcripts, $build_TE );
my $max_procs = 8;

( $build_index, $build_tRNAs, $build_rRNAs, $build_snRNAs, $build_miRNAs, $build_transcripts, $build_TE ) = (0,0,0,0,0,0,0);
( $min, $max, $mis, $misTE, $si_min, $si_max, $pi_min, $pi_max, $dir ) = ( 18, 29, 0, 3, 21, 21, 23, 29 );
$Pcheck ='true';

GetOptions (
	"fastq=s" => \@fastq,
	"fastq_n=s" => \@fastq_n,
	"dir=s" => \$dir,
	"min:i" => \$min,
	"max:i" => \$max,
	"si_min:i" => \$si_min,
	"si_max:i" => \$si_max,
	"pi_min:i" => \$pi_min,
	"pi_max:i" => \$pi_max,
	"mis:i" => \$mis,
	"misTE:i" => \$misTE,
	"html:s" => \$html_out,
	"PPPon:s" => \$Pcheck,
	"help"   =>  \$help,
	"ref:s" => \$ref,
	"tRNAs:s" => \$tRNAs,
	"rRNAs:s" => \$rRNAs,
	"snRNAs:s" => \$snRNAs,
	"miRNAs:s" => \$miRNAs,
	"transcripts:s" => \$transcripts,
	"TE:s" => \$TE,
	"build_index" => \$build_index,
	"build_tRNAs" => \$build_tRNAs,
	"build_snRNAs" => \$build_snRNAs,
	"build_miRNAs" => \$build_miRNAs,
	"build_transcripts" => \$build_transcripts,
	"build_rRNAs" => \$build_rRNAs,
	"build_TE" => \$build_TE
);

my $fq_collection = 'fastq_dir/';
mkdir $dir; mkdir $fq_collection;
$dir = $dir.'/' unless $dir =~ /\/$/;
mkdir $dir.'/css';mkdir $dir.'/js';
dircopy( $FindBin::Bin.'/css', $dir.'/css' );
dircopy( $FindBin::Bin.'/js', $dir.'/js' );

my $file = $dir.'report.txt';
open my $report, '>', $file or die "Cannot open $file $!\n";

my @toBuild = ( [$build_index, \$ref],  [$build_tRNAs, \$tRNAs], [$build_rRNAs, \$rRNAs], [$build_snRNAs, \$snRNAs], [$build_miRNAs, \$miRNAs], [$build_transcripts, \$transcripts], [$build_TE, \$TE] );
to_build ( \@toBuild, $report, $dir );

my $proc_child = ceil($max_procs / scalar(@fastq));
my $proc_grand_child = ceil($proc_child/4);
my $pm = Parallel::ForkManager->new($max_procs);
my $pm2 = Parallel::ForkManager->new($proc_grand_child);

$pm->run_on_finish( sub {
	my ($pid, $exit_code, $ident) = @_;
	print $report "Fastq fork $ident just finished ".
	"with PID $pid and exit code: $exit_code\n";
	die "Something went wrong!\n" if $exit_code != 0;
	});
$pm->run_on_start( sub {
	my ($pid,$ident)=@_;
	print $report "Fastq fork : $ident started, pid: $pid\n";
	});
$pm2->run_on_finish( sub {
	my ($pid, $exit_code, $ident) = @_;
	print $report "** Subgroup fork $ident just finished ".
	"with PID $pid and exit code: $exit_code\n";
	die "Something went wrong!\n" if $exit_code != 0;
	});
$pm2->run_on_start( sub {
	my ($pid,$ident)=@_;
	print $report "** Subgroup fork $ident started, pid: $pid\n";
	});


foreach my $child ( 0 .. $#fastq )
{
	my @suffix = ('.fastq', '.fastq.gz,', '.fq', '.fq.gz', 'ref', '.dat', '.fa','.fas','.fasta', '.txt');
	my ( $name, $path, $suffix ) = fileparse( $fastq[$child], @suffix );
	my ( $ref_name, $ref_path, $ref_suffix ) = fileparse( $ref, @suffix );
	my ( $TE_name, $TE_path, $TE_suffix ) = fileparse( $TE, @suffix );
	my ( $ex_name, $ex_path, $ex_suffix ) = fileparse( $transcripts, @suffix );

	$pm->start($fastq[$child]) and next;

	my $dir_fq = $dir.$fastq_n[$child].'/';
	mkdir $dir_fq;

	my $gen_dir = $dir_fq.'genome/';
	mkdir $gen_dir;

	my $size_dir = $dir_fq.'size/';
	mkdir $size_dir;

	my $fastq_resized = $dir_fq.$name.'_'.$min.'-'.$max.'.fastq';
	size_distribution (  $fastq[$child], $fastq_resized, $size_dir, $min, $max );

	my $sam_genome = $gen_dir.$fastq_n[$child].'_'.$min.'-'.$max.'.sam';
	my $sam_genome_unique = $gen_dir.$fastq_n[$child].'_'.$min.'-'.$max.'_unique.sam';
	my $fastq_prefix = $gen_dir.$fastq_n[$child].'_'.$min.'-'.$max;

	BWA_call ( $ref, $fastq_resized, $sam_genome, $mis, $proc_child, $report );
	my ( $fai_ref_hashP, $ma, $ma_uni ) = get_unique ( $sam_genome, $sam_genome_unique, $gen_dir, $fq_collection.$fastq_n[$child], 1, $report );

	die "No Reads mapped on the genome reference!\n" if $ma == 0;
	my $scale = 1000000 / $ma;
	sam_to_bam_bg ( $sam_genome_unique, $scale, $proc_child );
	sam_to_bam_bg ( $sam_genome, $scale, $proc_child );

	my $Gviz_dir = $gen_dir.'Gviz/';
	my $fai_file =  $gen_dir.'fai';
	mkdir $Gviz_dir;
	my $Gviz_dir_rand = $Gviz_dir.'rand/';
	mkdir $Gviz_dir_rand;
	my $Gviz_dir_uni = $Gviz_dir.'unique/';
	mkdir $Gviz_dir_uni;

	open my $gfai, '>', $fai_file;
	foreach my $k  ( sort keys %{$fai_ref_hashP} )
	{
		print $gfai "$k\t$fai_ref_hashP->{$k}\n";
	}
	close $gfai;
	bg_to_png ( $fai_file, $fastq_prefix.'_unique_plus.bedgraph', $fastq_prefix.'_unique_minus.bedgraph', $Gviz_dir_uni, 'Mb' );
	bg_to_png ( $fai_file, $fastq_prefix.'_plus.bedgraph', $fastq_prefix.'_minus.bedgraph', $Gviz_dir_rand, 'Mb' );

	my $group_dir = $dir_fq.'subgroups/';
	my $fastq_uni = $fq_collection.$fastq_n[$child].'_unique_mappers.fastq';
	my $fastq_all = $fq_collection.$fastq_n[$child].'_all_mappers.fastq';
	my ($bo, $mi, $pi) = subgroups ( $fastq_all, $group_dir, $mis, $misTE, $proc_child, $tRNAs, $rRNAs, $snRNAs, $miRNAs, $transcripts, $TE, $si_min, $si_max, $pi_min, $pi_max, $report);

	pie_chart($group_dir);

	open (my $dupnum, $gen_dir.'dup_mapnum.txt') || die "cannot open dup_mapnum.txt $!";
	my %dupnum_genome;
	my $header = <$dupnum>;
	while (<$dupnum>)
	{
		chomp $_;
		my @dupline = split /\t/, $_;
		$dupnum_genome{$dupline[0]} = [$dupline[1], $dupline[2]];
	}
	close $dupnum;

	my $mi_sam = $group_dir.'miRNAs.sam';
	mkdir $group_dir.'miRNAs/';
	my $mi_count_file =  $group_dir.'miRNAs/miRNAs_reads_counts.txt';
	my ( $mi_count, $mi_ref_size ) = sam_count ( $mi_sam );

	rpms_rpkm( $mi_count, $mi_ref_size, $ma, $mi_count_file, $pi, $mi, $bo );

	my (  $sam_transcripts, $sam_TEs ) = ( $group_dir.'transcripts.sam', $group_dir.'TEs.sam' );
	my @types = ($group_dir.'bonafide_reads.fastq', $group_dir.'miRNAs.fastq', $group_dir.'siRNAs.fastq', $group_dir.'piRNAs.fastq' );
	my @types_names = ('bonafide_reads', 'miRNAs', 'siRNAs', 'piRNAs');
	foreach my $grand_child ( 0 .. $#types )
	{
		my $type_dir = $group_dir.$types_names[$grand_child].'/';
		my $type_prefix = $types_names[$grand_child].'-';
		mkdir  $type_dir;
		$pm2->start($types[$grand_child]) and next;
		my ( $type_sam_genome, $type_sam_TEs, $type_sam_transcripts ) = ( $type_dir.$type_prefix.'genome.sam', $type_dir.$type_prefix.'TEs.sam', $type_dir.$type_prefix.'transcripts.sam' );
		my ( $type_sam_uni_genome, $type_sam_uni_TEs,  $type_sam_uni_transcripts ) = ( $type_dir.$type_prefix.'genome_unique.sam', $type_dir.$type_prefix.'TEs_unique.sam', $type_dir.$type_prefix.'transcripts_unique.sam' );
		my ( $type_uni_genome_fastq, $type_uni_TEs_fastq,  $type_uni_transcripts_fastq ) = ( $fq_collection.$fastq_n[$child].'-'.$type_prefix.'genome_uni.fastq', $fq_collection.$fastq_n[$child].'-'.$type_prefix.'TEs_uni.fastq', $fq_collection.$fastq_n[$child].'-'.$type_prefix.'transcripts_uni.fastq');
		my ( $type_genome_fastq, $type_TEs_fastq,  $type_transcripts_fastq ) = ( $fq_collection.$fastq_n[$child].'-'.$type_prefix.'genome.fastq', $fq_collection.$fastq_n[$child].'-'.$type_prefix.'TEs.fastq', $fq_collection.$fastq_n[$child].'-'.$type_prefix.'transcripts.fastq');
		my $type_sequence_hashP = get_fastq_seq ( $types[$grand_child] );

		if ( $grand_child == 1 )
		{
			BWA_call ( $TE, $types[$grand_child],  $type_sam_TEs, $misTE, $proc_child, $report );
			BWA_call ( $transcripts, $types[$grand_child], $type_sam_transcripts, $mis, $proc_child, $report );
		 	BWA_call ( $ref, $types[$grand_child], $type_sam_genome, $mis, $proc_child, $report );
			extract_sam ( undef, $type_sam_TEs, $type_sam_TEs, $type_sam_uni_TEs, $type_uni_TEs_fastq, $type_uni_TEs_fastq );
			extract_sam ( undef, $type_sam_transcripts, $type_sam_transcripts, $type_sam_uni_transcripts, $type_transcripts_fastq, $type_uni_transcripts_fastq );
			extract_sam ( undef, $type_sam_genome, $type_sam_genome, $type_sam_uni_genome, $type_genome_fastq, $type_uni_genome_fastq );
		}
		else
		{
			extract_sam ( $type_sequence_hashP, $sam_TEs, $type_sam_TEs, $type_sam_uni_TEs, $type_TEs_fastq, $type_uni_TEs_fastq );
			extract_sam ( $type_sequence_hashP, $sam_transcripts, $type_sam_transcripts, $type_sam_uni_transcripts, $type_transcripts_fastq, $type_uni_transcripts_fastq );
			extract_sam ( $type_sequence_hashP, $sam_genome, $type_sam_genome, $type_sam_uni_genome, $type_genome_fastq, $type_uni_genome_fastq );
		}

		my $ex_count_file =  $type_dir.$type_prefix.'transcripts_reads_counts.txt';
		my ( $ex_count, $ex_ref_size ) =  sam_count_mis ( $type_sam_transcripts );
		rpms_rpkm_te( $ex_count, $ex_ref_size, $ma, $ex_count_file, $pi, $mi, $bo );

		my ( $TEs_count, $TEs_ref_size, $TEs_count_NoM, $TEs_count_M ) = sam_count_mis ( $type_sam_TEs );
		my $TEs_count_file = $type_dir.$type_prefix.'TEs_reads_counts.txt';
		my $TEs_count_file_M = $type_dir.$type_prefix.'TEs_reads_counts_mismatches.txt';
		my $TEs_count_file_noM = $type_dir.$type_prefix.'TEs_reads_counts_nomismatches.txt';
		rpms_rpkm_te( $TEs_count, $TEs_ref_size, $ma, $TEs_count_file, $pi, $mi, $bo );
		rpms_rpkm_te( $TEs_count_NoM, $TEs_ref_size, $ma, $TEs_count_file_noM, $pi, $mi, $bo );
		rpms_rpkm_te( $TEs_count_M, $TEs_ref_size, $ma, $TEs_count_file_M, $pi, $mi, $bo );

		sam_to_bam_bg ( $type_sam_TEs, $scale, $grand_child );
		sam_sorted_bam ( $type_sam_transcripts, $grand_child ); sam_sorted_bam ( $type_sam_uni_transcripts, $grand_child ); 
		sam_sorted_bam ( $type_sam_uni_TEs, $grand_child ); 

		my $Gviz_TEs =  $type_dir.'Gviz_TEs/';
		mkdir $Gviz_TEs;
		bg_to_png ( $group_dir.'TEs.fai', $type_dir.$type_prefix.'TEs_plus.bedgraph', $type_dir.$type_prefix.'TEs_minus.bedgraph', $Gviz_TEs, 'Kb' );

		my $Gviz_genome=  $type_dir.'Gviz_genome/';
		my $Gviz_genome_rand = $Gviz_genome.'rand/';
		my $Gviz_genome_uni = $Gviz_genome.'unique/';
		mkdir $Gviz_genome; mkdir $Gviz_genome_uni; mkdir $Gviz_genome_rand;

		sam_to_bam_bg ( $type_sam_genome, $scale, $grand_child );
		sam_to_bam_bg ( $type_sam_uni_genome, $scale, $grand_child );

		bg_to_png ( $fai_file, $type_dir.$type_prefix.'genome_unique_plus.bedgraph', $type_dir.$type_prefix.'genome_unique_minus.bedgraph', $Gviz_genome_uni, 'Mb' );
		bg_to_png ( $fai_file, $type_dir.$type_prefix.'genome_plus.bedgraph', $type_dir.$type_prefix.'genome_minus.bedgraph', $Gviz_genome_rand, 'Mb' );

		#HTML Details
		my $prefix_details_pages = $dir.$fastq_n[$child].'-'.$types_names[$grand_child];
		details_pages ( $type_dir, $prefix_details_pages, \@fastq_n, $fastq_n[$child], $misTE, $dir, $Pcheck );

		$pm2->finish();
	}
	$pm2->wait_all_children;

	if ( $Pcheck eq 'true' )
	{
		my $ppp = $group_dir.'PPPartners/'; mkdir $ppp;
		print $report "ping_pong_partners $group_dir/piRNAs/TEs.sam $ppp\n";
		ping_pong_partners ( $group_dir.'TEs.fai', $group_dir.'piRNAs/piRNAs-TEs_sorted.bam', $ppp, $pi_min );
		my $ppp_page = $dir.$fastq_n[$child].'-piRNAs-PPP.html';
		ppp_page ( $group_dir, $ppp_page, \@fastq_n, $fastq_n[$child], $ppp, $dir );
	}

	#HTML Main Webpage
	my $index_page = $dir.$fastq_n[$child].'.html';
	main_page ( $gen_dir, $index_page, \@fastq_n, $fastq_n[$child], $ma, $ma_uni, $dir );
	copy ($index_page, $html_out) if $child == 0;
	#HTML Menu
	my $menu_page = $dir.$fastq_n[$child].'-sub.html';
	menu_page ( $group_dir, $menu_page, \@fastq_n, $fastq_n[$child], $min, $max, $si_min, $si_max, $pi_min, $pi_max, $dir );
	unlink glob "'$group_dir'*.sam"; unlink glob "'$group_dir'*.fastq";
	$pm->finish(); # pass an exit code to finish
}
$pm->wait_all_children;
unlink glob "'$dir'"."dataset_*symlink.fa*";
print $report "Job done!\n";
close $report;

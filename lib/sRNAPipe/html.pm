package sRNAPipe::html;

use strict;
use warnings;
use File::Basename;
use File::Copy::Recursive qw( dircopy );

use Exporter;
our @ISA = qw( Exporter );
our @EXPORT_OK = qw( &main_page &details_pages &menu_page &ppp_page &copy_css &copy_js );

sub main_page
{
  my ( $dir, $file, $list_mainTabP, $current, $ma, $ma_uni, $dir_root ) = @_;
  my ( $futHashP, $uniqueTabP, $randTabP, $pngTabP ) = get_genome ( $dir, $dir_root );

  open my $h, '>', $file || die "cannot create $file $!\n";
  header ( $h );
  navbar ( $h, $list_mainTabP, $current );
  print $h "<div class=\"container\"><p><a class=\"btn\" href=\"$current-sub.html\">View details &raquo;</a></p></div>\n";
  futurette( $h, $current, $pngTabP, $futHashP );
  print  $h "<div class=\"container\"><h2>mappers #: $ma</h2><h2>unique mappers #: $ma_uni</h2> </div>\n";
  carousel2( $h, $uniqueTabP, $randTabP, $dir_root );
  footer($h);
  close $h;
}

sub menu_page
{
  my ( $dir, $file, $list_mainTabP, $current, $min, $max, $simin, $simax, $pimin, $pimax, $dir_root ) = @_;
  my $html_ref = $1 if $dir =~ /$dir_root(.*)/;
  open my $h, '>', $file || die "cannot create $file $!\n";
  header($h);
  navbar ( $h, $list_mainTabP, $current );
  span( $h, $current, $min, $max, $simin, $simax, $pimin, $pimax );
  print $h "  <div class=\"container\"> <div class=\"row text-center\">  <img  src=\"$html_ref/pie_chart.png\"/><br />\n";
  print $h "  <A HREF=\"$html_ref/repartition.txt\">text file</A><br/>\n </div></div>";
  footer($h);
  close $h;
}

sub details_pages
{
  my ( $dir_details, $prefix, $list_mainTabP, $current, $misTE, $dir_root, $ppp ) = @_;
  my ($Hex, $HTE, $HG, $NonUniTE, $NonUniG, $UniG ) = get_subgroups( $dir_details, $current, $misTE, $dir_root );

  my $html_ref = $1.'-PPP.html' if $prefix =~ /$dir_root(.*)/;
  open my $h, '>',  $prefix.'-TEs.html' || die "cannot create  $prefix-TEs.html $!\n";
  header($h);
  navbar ( $h, $list_mainTabP, $current );
  if ( $prefix =~ /piRNAs$/ && $ppp eq 'true' )
  {
    print $h " <div class=\"container\">";
    print $h " <p><a class=\"btn\" href=\"$html_ref\">Ping Pong Partners</a></p>\n";
    print $h "</div>";
  }
  fut($h,'Transposable elements',$HTE);
  carousel($h,$NonUniTE,$dir_root);
  footer($h);
  close $h;

  open $h, '>',  $prefix.'-genome.html' || die "cannot create  $prefix-genome.html $!\n";
  header($h);
  navbar ( $h, $list_mainTabP, $current );
  fut($h,'Genome',$HG);
  carousel2($h,$UniG, $NonUniG,$dir_root);
  footer($h);
  close $h;

  open  $h, '>',  $prefix.'-transcripts.html' || die "cannot create  $prefix-transcripts.html $!\n";
  header($h);
  navbar ( $h, $list_mainTabP, $current );
  fut($h,'transcripts',$Hex);
  footer($h);
  close $h;
}

sub ppp_page
{
  my ( $dir, $file, $list_mainTabP, $current, $ppp, $dir_root ) = @_;

  my $ppp_file = $ppp.'ppp.txt';
  open my $h, '>', $file || die "cannot create $file $!\n";
  header($h);
  navbar ( $h, $list_mainTabP, $current );
  print $h '<div class="container"> <table class="wb-tables table table-striped table-hover">'."\n";
  print $h '<thead>
  <tr>
    <th data-sortable="true">ID</th>
    <th data-sortable="true">overlap sum</th>
    <th data-sortable="true">ten overlap sum</th>
    <th data-sortable="true">mean</th>
    <th data-sortable="true">standard deviation</th>
  	<th data-sortable="true">z-score</th>
  	<th data-sortable="true">p-value</th>
	</tr>
  </thead>
  <tbody>';

  open my $f, '<', $ppp_file || die "cannot open $ppp_file  $!\n";
  while ( <$f> )
  {
    chomp;
    print $h '<tr>';
    my ( $id, $sum, $ten, $mean, $sd, $zscore, $prob) = split /\t/, $_;
    if( -d "$ppp/$id" )
    {
      my $sub_html = $ppp.$id.'.html';
      my $sub_html_ref = $1.$id if $ppp =~ /$dir_root(.*)/;
      print $h "<td> <a href=\"$sub_html_ref.html\">$id</a> </td>";

      open my $sub, '>', $sub_html || die "cannot create $sub_html\n";
      {
        header($sub);
        print $sub "
					<div align=\"center\">
					<h2>$id</h2>
					<p> <img class=\"featurette-image\" src=\"$id/histogram.png\"/></p>
					<p><a href=\"$id/overlap_size.txt\">ping pong signature</a></p>
					<p><a href=\"$id/sensPPP.txt\">sense reads with PPP</a></p>
					<p><a href=\"$id/antisensPPP.txt\">reverse reads  with PPP</a></p>
					<p><a href=\"$id/sens.txt\">sense reads without PPP</a></p>
					<p><a href=\"$id/antisens.txt\">reverse reads without PPP</a></p>
					</div>";
        footer($sub);
      }
      close $sub;
    }
    else { print $h "<td> $id </td>\n"; }
    print $h "<td> $sum </td><td> $ten </td><td> $mean </td><td> $sd </td><td> $zscore </td><td> $prob </td>\n";

    print $h '</tr>';
  }
  close $f;
  print $h "</tbody></table></div>";
  footer($h);
  close $h;
}

sub get_genome
{
  my ( $dir, $dir_root ) = @_;
  my ( %hash, @group,  @Unique, @NonUnique, @png );

  my $fut = "'$dir'".'/*';
  my @fut = glob $fut;


  foreach my $fr ( @fut )
  {
    my $f = $1 if $fr =~ /$dir_root(.*)/;
    if ( $fr =~ /.*Gviz/ )
    {
      my $nu = "'$fr'".'/rand/*';
      @NonUnique =  glob $nu;
      my $u = "'$fr'".'/unique/*';
      @Unique =  glob $u;
    }
    elsif ( $f =~ /.*distribution\.txt$/ ) { $hash{'mappers size distribution (txt)'} = $f; }
    elsif ( $f =~ /.*distribution\.png$/ ) { push @png, $f; }
    elsif ( $f =~ /.*unique\.fastq$/ ) { $hash{'unique mappers (fastq.gz)'} = $f.'.gz'; `gzip '$fr'`; }
    elsif ( $f =~ /.*rejected\.fastq$/ ) { $hash{'unmapped (fastq.gz)'} = $f.'.gz'; `gzip '$fr'`; }
    elsif ( $f =~ /.*all\.fastq$/ ) { $hash{'mappers (fastq.gz)'} = $f.'.gz'; `gzip '$fr'`; }
    elsif ( $f =~ /.*dup_unique\.txt$/ ) { $hash{'unique mappers (txt)'} = $f; }
    elsif ( $f =~ /.*dup_mapnum\.txt$/ ) { $hash{'mappers (txt)'} = $f; }
    elsif ( $f =~ /.*dup_nonmapp\.txt$/ ) { $hash{'unmapped (txt)'} = $f; }
    elsif ( $f =~ /.*_unique_sorted\.bam$/ ) { $hash{'unique alignment (bam)'} = $f; }
    elsif ( $f =~ /.*_sorted\.bam$/ ) { $hash{'alignment (bam)'} = $f; }
    elsif ( $f =~ /.*unique_plus.bedgraph/) { $hash{'bedgraph unique plus strand'} = $f; }
    elsif ( $f =~ /.*unique_minus.bedgraph/) { $hash{'bedgraph unique minus strand'} = $f; }
    elsif ( $f =~ /.*plus.bedgraph/) { $hash{'bedgraph plus strand'} = $f; }
    elsif ( $f =~ /.*minus.bedgraph/) { $hash{'bedgraph minus strand'} = $f; }
    else { unlink $fr; }
  }
  return (\%hash, \@Unique, \@NonUnique, \@png);
}

sub span
{
  my ( $file, $name, $min, $max, $simin, $simax, $pimin, $pimax ) = @_;

  print $file "
<div class=\"container  text-center\">
  <div class=\"row-fluid\">
      <div class=\"span6\">
        <h2>Bonafide</h2>
        reads of size between $min and $max<br>with no mi, sn, t and r RNAs
        <p><a class=\"btn\" href=\"$name-bonafide_reads-genome.html\">Genome</a></p>
        <p><a class=\"btn\" href=\"$name-bonafide_reads-TEs.html\">TE</a></p>
        <p><a class=\"btn\" href=\"$name-bonafide_reads-transcripts.html\">Transcripts</a></p>
        <div class=\"row-fluid\">
          <div class=\"span6\">
            <h2>siRNAs</h2>
            bonafide reads of size between $simin and $simax
            <p><a class=\"btn\" href=\"$name-siRNAs-genome.html\">Genome</a></p>
            <p><a class=\"btn\" href=\"$name-siRNAs-TEs.html\">TE</a></p>
            <p><a class=\"btn\" href=\"$name-siRNAs-transcripts.html\">Transcripts</a></p>
          </div>
          <div class=\"span6\">
            <h2>piRNAs</h2>
            bonafide reads of size between $pimin and $pimax
            <p><a class=\"btn\" href=\"$name-piRNAs-genome.html\">Genome</a></p>
            <p><a class=\"btn\" href=\"$name-piRNAs-TEs.html\">TE</a></p>
            <p><a class=\"btn\" href=\"$name-piRNAs-transcripts.html\">Transcripts</a></p>
          </div>
        </div>
      </div>
    <div class=\"span6\">
       <h2>miRNAs</h2>
       <p><a class=\"btn\" href=\"$name-miRNAs-genome.html\">Genome</a></p>
       <p><a class=\"btn\" href=\"$name-miRNAs-TEs.html\">TE</a></p>
       <p><a class=\"btn\" href=\"$name-miRNAs-transcripts.html\">Transcripts</a></p>
    </div>
  </div>
</div>
";
}

sub get_subgroups
{
  my ( $dir, $name, $misTE, $dir_root ) = @_;
  my (%Hex, %HTE, %HG, @group, @png, @pngTE,  @NonUniTE, @UniG, @NonUniG );

  my $fut = "'$dir'".'/*';
  my @fut = glob $fut;
  my $f ='';
  foreach my $fr ( @fut )
  {
    $f = $1 if $fr =~  /$dir_root(.*)/;

    if ( $f =~ /genome_unique_sorted\.bam$/ ) { $HG{'genome unique mappers (sorted bam)'} =  $f; }
    elsif ( $f =~ /genome_sorted\.bam$/ ) { $HG{'genome mappers (sorted bam)'} = $f; }
    elsif ( $f =~ /miRNAs_reads_counts\.txt$/ ) { $HG{'miRNAs per type (txt)'} = $f; }
    elsif ( $f =~ /genome_unique_plus\.bedgraph$/) { $HG{'bedgraph unique plus strand'} = $f; }
    elsif ( $f =~ /genome_unique_minus\.bedgraph$/) { $HG{'bedgraph unique minus strand'} = $f; }
    elsif ( $f =~ /genome_plus\.bedgraph$/) { $HG{'bedgraph plus strand'} = $f; }
    elsif ( $f =~ /genome_minus\.bedgraph$/) { $HG{'bedgraph minus strand'} = $f; }
    elsif ( $f =~ /TEs_plus\.bedgraph$/) { $HTE{'bedgraph plus strand'} = $f; }
    elsif ( $f =~ /TEs_minus\.bedgraph$/) { $HTE{'bedgraph minus strand'} = $f; }
    elsif ( $f =~ /transcripts_sorted\.bam$/) { $Hex{'transcripts mappers (sorted bam)'} = $f;}
    elsif ( $f =~ /transcripts_unique_sorted\.bam$/) { $Hex{'transcripts unique mappers (sorted bam)'} = $f;}
    elsif ( $f =~ /transcripts_reads_counts\.txt$/) { $Hex{'read number per transcript (txt)'} = $f;}
    elsif ( $f =~ /TEs_reads_counts\.txt$/) { $HTE{"read number per TE 0 to $misTE mismatches (txt)"} = $f; }
    elsif ( $f =~ /TEs_reads_counts_mismatches\.txt$/) { $HTE{"read number per TE with 1 to $misTE mismatches (txt)"} = $f; }
    elsif ( $f =~ /TEs_reads_counts_nomismatches\.txt$/) { $HTE{'read number per TE with no mismatch (txt)'} = $f; }
    elsif ( $f =~ /TEs_unique_sorted\.bam$/) { $HTE{'TEs unique mappers (sorted bam)'} = $f; }
    elsif ( $f =~ /TEs_sorted\.bam$/) { $HTE{'TEs mappers (sorted bam)'} = $f; }
    elsif ( $fr =~ /.*Gviz_TEs/ )
    {
      my $nu = "'$fr'".'/*';
      @NonUniTE =  glob $nu;
    }
    elsif ( $fr =~ /.*Gviz_genome/ )
    {
      my $nu = "'$fr'".'/rand/*';
      @NonUniG =  glob $nu;
      my $u = "'$fr'".'/unique/*';
      @UniG =  glob $u;
    }
    else { unlink $fr; }
  }
  return (\%Hex, \%HTE, \%HG, \@NonUniTE,  \@NonUniG, \@UniG);
}

sub header
{
  my $file = shift;
  print $file "
  <!DOCTYPE html>
  <html lang=\"en\">
  <head>
  <meta charset=\"utf-8\">
  <title>pipeline</title>
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
  <meta name=\"description\" content=\"\">
  <meta name=\"author\" content=\"\">
  <!-- Le styles -->
  <link href=\"css/bootstrap.css\" rel=\"stylesheet\">
  <link href=\"css/bootstrap-table.css\" rel=\"stylesheet\">
  <style type=\"text/css\">
  body {
    padding-top: 60px;
    padding-bottom: 40px;
  }
  div#page {
    width: 940px;
    background-color: #fff;
    margin: 0 auto;
    text-align: left;
    border-color: #fff;
    border-style: none solid solid;
    border-width: medium 1px 1px;
  }

  div.content {
   	display: none;
	  float: right;
	  width: 550px;
  }
  div.content a, div.navigation a {
    text-decoration: none;
   	color: #777;
  }
  div.content a:focus, div.content a:hover, div.content a:active {
    text-decoration: underline;
  }
  
  div.controls {
    margin-top: 5px;
	  height: 23px;
  }
  
  div.controls a {
	  padding: 5px;
  }
  div.ss-controls {
	  float: left;
  }
  div.nav-controls {
	  float: right;
  }
  div.slideshow-container {
	  position: relative;
	  clear: both;
	  height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
  }
  div.loader {
  	position: absolute;
  	top: 0;
  	left: 0;
    background-image: url('loader.gif');
    background-repeat: no-repeat;
    background-position: center;
  	width: 550px;
  	height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
  }
  div.slideshow {
    
  }
  
  div.slideshow span.image-wrapper {
	  display: block;
	  position: absolute;
	  top: 0;
	  left: 0;
  }
  div.slideshow a.advance-link {
	  display: block;
 	  width: 550px;
    height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
    line-height: 502px; /* This should be set to be at least the height of the largest image in the slideshow */
    text-align: center;
  }
  div.slideshow a.advance-link:hover, div.slideshow a.advance-link:active, div.slideshow a.advance-link:visited {
    text-decoration: none;
  }
  div.slideshow img {
    vertical-align: middle;
	border: 1px solid #ccc;
  }

  div.image-title {
    font-weight: bold;
    font-size: 1.4em;
  }

  div.image-desc {
    line-height: 1.3em;
    padding-top: 12px;
  }
  div.navigation {

  }
  ul.thumbs {
	clear: both;
	margin: 0;
	padding: 0;
  }
  ul.thumbs li {
	float: none;
	padding: 0;
	margin: 0;
    list-style: none;
  }
  a.thumb {
	padding: 0;
	display: inline;
	border: none;
  }
  ul.thumbs li.selected a.thumb {
	color: #000;
    font-weight: bold;
  }
  a.thumb:focus {
	outline: none;
  }
  ul.thumbs img {
	border: none;
	display: block;
  }
  div.pagination {
	clear: both;
  }
  div.navigation div.top {
    margin-bottom: 12px;
	height: 11px;
  }
  div.navigation div.bottom {
    margin-top: 12px;
  }
  div.pagination a, div.pagination span.current, div.pagination span.ellipsis {
	display: block;
	float: left;
    margin-right: 2px;
	padding: 4px 7px 2px 7px;
	border: 1px solid #ccc;
  }
  div.pagination a:hover {
    background-color: #eee;
    text-decoration: none;
  }
  div.pagination span.current {
    font-weight: bold;
    background-color: #000;
    border-color: #000;
	color: #fff;
  }
  div.pagination span.ellipsis {
	border: none;
	padding: 5px 0 3px 2px;
  }
  
  div.download {
	float: right;
  }
  
  div.caption-container {
	position: relative;
	clear: left;
	height: 75px;
  }
  span.image-caption {
	display: block;
	position: absolute;
	width: 550px;
	top: 0;
	left: 0;
  }
  div.caption {
	padding: 12px;
  }

  /* Featurettes
  ------------------------- */

  .featurette {
  padding-top: 20px; /* Vertically center images part 1: add padding above and below text. */
  overflow: hidden; /* Vertically center images part 2: clear their floats. */
  text-align: center;
  }

  .featurette-p
  {
   text-align: left;
  }

  .featurette-image {
  margin-top: 10px; /* Vertically center images part 3: negative margin up the image the same amount of the padding to center it. */
  width: 450px;
  height: auto;
  }

  </style>
  <link href=\"css/bootstrap-responsive.css\" rel=\"stylesheet\">
  </head>
  <body>
  ";
}

sub navbar
{
  my ( $file, $fastq, $actif ) = @_;

  print $file "
  <div class=\"navbar navbar-inverse navbar-fixed-top\">
  <div class=\"navbar-inner\">
  <div class=\"container\">
  <button type=\"button\" class=\"btn btn-navbar\" data-toggle=\"collapse\" data-target=\".nav-collapse\">
  <span class=\"icon-bar\"></span>
  <span class=\"icon-bar\"></span>
  <span class=\"icon-bar\"></span>
  </button>
  <a class=\"brand\" href=\"report.txt\">Report</a>
  <div class=\"nav-collapse collapse\">
  <ul class=\"nav\">
  ";
  for (my $i = 0 ; $i <= $#{$fastq}; $i++)
  {
    # my $fa = basename($fastq->[$i],'.dat');
    my $fa = $fastq->[$i];
    if ($actif eq $fa){  print $file "<li class=\"active\"><a href=\"$fastq->[$i].html\">$fa</a></li>";}
    else {print $file "<li><a href=\"$fastq->[$i].html\">$fa</a></li>" ;}
  }
  print $file "
  </ul>
  </div><!--/.nav-collapse -->
  </div>
  </div>
  </div>";
}

sub footer
{
  my $file = shift;
  print $file "
  <!-- FOOTER -->
  <div class=\"container\">
  <footer>
  
  </footer>
  </div>
  <!-- Le javascript
  ================================================== -->
  <!-- Placed at the end of the document so the pages load faster -->
  <script type=\"text/javascript\" src=\"js/filter.js\"></script>
  <script type=\"text/javascript\" src=\"js/jquery.js\"></script>
  <script type=\"text/javascript\" src=\"js/jquery-1.3.2.js\"></script>
  <script type=\"text/javascript\" src=\"js/jquery.galleriffic.js\"></script>
  <script type=\"text/javascript\" src=\"js/jquery.opacityrollover.js\"></script>
  <script type=\"text/javascript\" src=\"js/bootstrap-table.js\"></script>
  <script type=\"text/javascript\" src=\"js/bootstrap.min.js\"></script>
  <script type=\"text/javascript\">
  jQuery(document).ready(function(\$) {
    // We only want these styles applied when javascript is enabled
    \$('div.navigation').css({'width' : '300px', 'float' : 'left'});
    \$('div.content').css('display', 'block');
    
    \$(\".each-gallery\").each(function(i){
      // Initially set opacity on thumbs and add
      // additional styling for hover effect on thumbs
        var onMouseOutOpacity = 0.67;
      \$('#thumbs + i + ul.thumbs li').opacityrollover({
      mouseOutOpacity:   onMouseOutOpacity,
      mouseOverOpacity:  1.0,
      fadeSpeed:         'fast',
      exemptionSelector: '.selected'
      });
      
      // Initialize Advanced Galleriffic Gallery
      var gallery = \$('#thumbs'+i).galleriffic({
      delay:                     2500,
      numThumbs:                 22,
      preloadAhead:              10,
      enableTopPager:            true,
      enableBottomPager:         true,
      maxPagesToShow:            7,
      imageContainerSel:         '#slideshow'+ i,
      controlsContainerSel:      '#controls' + i,
      captionContainerSel:       '#caption' + i,
      loadingContainerSel:       '#loading' + i,
      renderSSControls:          true,
      renderNavControls:         true,
      playLinkText:              'Play',
      pauseLinkText:             'Pause',
      prevLinkText:              '&lsaquo; Previous',
      nextLinkText:              'Next &rsaquo;',
      nextPageLinkText:          'Next &rsaquo;',
      prevPageLinkText:          '&lsaquo; Prev',
      enableHistory:             false,
      autoStart:                 false,
      syncTransitions:           true,
      defaultTransitionDuration: 900,
      onSlideChange:             function(prevIndex, nextIndex) {
        // 'this' refers to the gallery, which is an extension of \$('#thumbs')
        this.find('ul.thumbs').children()
        .eq(prevIndex).fadeTo('fast', onMouseOutOpacity).end()
        .eq(nextIndex).fadeTo('fast', 1.0);
      },
      onPageTransitionOut:       function(callback) {
        this.fadeTo('fast', 0.0, callback);
      },
      onPageTransitionIn:        function() {
        this.fadeTo('fast', 1.0);
      }
      });
    });
  });
  </script>
  </body>
  </html>
  ";
}

sub carousel
{
  my ($file, $non_unique, $dir_root) = @_;
  my $ac = 0;
  print $file "
  <div id=\"page\">
  <div id=\"container\">
  <div class=\"each-gallery\">
  <div id=\"gallery\" class=\"content\">
  <div id=\"controls0\" class=\"controls\"></div>
  <div class=\"slideshow-container\">
  <div id=\"loading0\" class=\"loader\"></div>
  <div id=\"slideshow0\" class=\"slideshow\"></div>
  </div>
  <div id=\"caption0\" class=\"caption-container\">Reads randomly assigned</div>
  </div>
  <div id=\"thumbs0\" class=\"navigation\">
  <input type=\"text\" id=\"myInput0\" onkeyup=\"search(this)\" placeholder=\"Search for names...\">
  <ul class=\"thumbs noscript\">
  ";
  foreach my $u (@{$non_unique})
  {
    my $name = basename($u,'.png');
    $u = $1 if $u =~ /$dir_root(.*)/;
    print $file "
    <li>
    <a class=\"thumb\"  href=\"$u\" title=\"$name\">$name</a>
    </li>
  	";
  }
  print $file "
  </ul>
  </div>
  <div style=\"clear: both;\"></div></div>
  </div>
  </div>
  ";
}

sub carousel2
{
  my ($file, $unique, $non_unique, $dir_root) = @_;
  print $file "
  <div id=\"page\">
  <div id=\"container\">
  <div class=\"each-gallery\">
  <div id=\"gallery\" class=\"content\">
  <div id=\"controls0\" class=\"controls\"></div>
  <div class=\"slideshow-container\">
  <div id=\"loading0\" class=\"loader\"></div>
  <div id=\"slideshow0\" class=\"slideshow\"></div>
  </div>
  <div id=\"caption0\" class=\"caption-container\">Uniquely mapped reads</div>
  </div>
  <div id=\"thumbs0\" class=\"navigation\">
  <input type=\"text\" id=\"myInput0\" onkeyup=\"search(this)\" placeholder=\"Search for names...\">
  <ul class=\"thumbs noscript\">
  ";

  foreach my $u (@{$unique})
  {
    my $name = basename($u,'.png');
    $u = $1 if $u =~ /$dir_root(.*)/;
    print $file "
    <li>
    <a class=\"thumb\"  href=\"$u\" title=\"$name\">$name</a>
    </li>
    ";
  }
  print $file "
  </ul>
  </div>
  </div>
  <div id=\"page\">
  <div id=\"container\">
  <div class=\"each-gallery\">
  <div id=\"gallery\" class=\"content\">
  <div id=\"controls1\" class=\"controls\"></div>
  <div class=\"slideshow-container\">
  <div id=\"loading1\" class=\"loader\"></div>
  <div id=\"slideshow1\" class=\"slideshow\"></div>
  </div>
  <div id=\"caption1\" class=\"caption-container\">Reads randomly assigned</div>
  </div>
  <div id=\"thumbs1\" class=\"navigation\">
  <input type=\"text\" id=\"myInput1\" onkeyup=\"search(this)\" placeholder=\"Search for names...\">
  <ul class=\"thumbs noscript\">
  ";

  foreach my $nu (@{$non_unique})
  {
    my $name = basename($nu,'.png');
    $nu = $1 if $nu =~ /$dir_root(.*)/;
    print $file "
        <li>
        <a class=\"thumb\"  href=\"$nu\" title=\"$name\">$name</a>
        </li>
        ";
  }
  print $file "
  </ul>
  </div>
  <div style=\"clear: both;\"></div></div>
  </div>
  </div>
  ";
}

sub futurette
{
  my ($file, $name, $png, $hash) = @_;
  print $file "
  <div class=\"container\">
  <div class=\"featurette\">
  <h1>$name</h1>
  <p class=\"featurette-p\">
  ";
  foreach my $k (sort  keys %{$hash})
  {
    print $file "<A HREF=\"".${$hash}{$k}."\">$k</A><br/> \n" ;
  }

  print $file "
  </p>";

  foreach my $pn (@{$png}){print $file "<img  class=\"featurette-image\"  src=\"$pn\"/><br />";}

  print $file "
  </div>
  </div>
  ";
}

sub fut
{
  my ($file, $name, $hash) = @_;
  print $file "
  <div class=\"container\">
  <div class=\"featurette\">
  <h1>$name</h1>
  <p class=\"featurette-p\">
  ";

  foreach my $k (sort { ${$hash}{$a} cmp ${$hash}{$b} } keys %{$hash})
  {
    print $file "<A HREF=\"".${$hash}{$k}."\">$k</A><br/> \n" ;
  }

  print $file "
  </p>
  </div>
  </div>
  ";
}

sub get_distri_exon
{
  my ($dir, $name) = @_;
  my (@out,@group);
  my $group = "'$dir'".'/'."'$name'".'-subgroups-bonafide_reads-transcripts-*distribution-*.png';
  @group = glob $group;
  foreach (my $g =0; $g <= $#group; $g++)
  {
    if ($group[$g] =~ /.*($name-subgroups-bonafide_reads-transcripts-.*distribution-.*\.png)/ )
    {
      my $tmp = $1;
      push @out, $1;
    }
  }
  return (\@out);
}

sub get_distri_TE
{
  my ($dir, $name) = @_;
  my (@out,@group);
  my $group = "'$dir'".'/'."'$name'".'-subgroups-bonafide_reads-TE-*distribution-*.png';
  @group = glob $group;
  foreach (my $g =0; $g <= $#group; $g++)
  {
    if ($group[$g] =~ /.*($name-subgroups-bonafide_reads-TE-.*distribution-.*\.png)/ )
    {
      my $tmp = $1;
      push @out, $1;
    }
  }
  return (\@out);
}

sub get_PPP
{
  my ($dir,$name) = @_;
  my (%distri,@group);
  my $group = "'$dir'".'/'."'$name'".'-subgroups-bonafide_reads-TE-PPPartners-*';
  @group = glob $group;

  foreach (my $g =0; $g <= $#group; $g++)
  {
    if ($group[$g] =~ /.*($name-subgroups-bonafide_reads-TE-PPPartners-.*)/ )
    {
      my $tmp = $1;
      if ($tmp =~ /PPPartners-(.*?)-sens\.txt$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[0] = $tmp;
      }
      elsif ($tmp =~ /PPPartners-(.*?)-antisens\.txt$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[1] = $tmp;
      }
      elsif ($tmp =~ /PPPartners-(.*?)-sensPPP\.txt$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[2] = $tmp;
      }
      elsif ($tmp =~ /PPPartners-(.*?)-antisensPPP\.txt$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[3] = $tmp;
      }
      elsif ($tmp =~ /PPPartners-(.*?)-overlap_size\.txt$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[4] = $tmp;
      }
      elsif ($tmp =~ /PPPartners-(.*?)-histogram\.png$/)
      {
        $distri{$1} = ['','','','','',''] unless exists $distri{$1};
        $distri{$1}->[5] = $tmp;
      }
    }
  }
  return \%distri;
}

sub PPPrint
{
  my ($h, $hash) = @_;
  my $cmp = 0;

  print $h "<div class=\"container\">\n";
  print $h "<div class=\"row text-center\">";
  while ( my ($k,$v) = each %{$hash} )
  {
    print $h "</div><div class=\"row text-center\">" if $cmp != 0 && $cmp % 2 == 0;
    print $h "
    
    <div class=\"span6\">
    <h2>$k</h2>
    <p class=\"featurette-p\"> <img src=\"$v->[5]\"/></p>
    <p class=\"featurette-p\"><a href=\"$v->[4]\">ping pong signature</a></p>
    <p class=\"featurette-p\"><a href=\"$v->[2]\">sense reads with PPP</a></p>
    <p class=\"featurette-p\"><a href=\"$v->[3]\">reverse reads  with PPP</a></p>
    <p class=\"featurette-p\"><a href=\"$v->[0]\">sense reads without PPP</a></p>
    <p class=\"featurette-p\"><a href=\"$v->[1]\">reverse reads without PPP</a></p>
    </div>
    ";
    $cmp++;
  }

  print $h "</div></div>";
}

sub printDistri
{
  my ($h, $tab) = @_;
  my ($txt, $name);
  my $cmp = 0;
  print $h "<div class=\"container\">\n";
  print $h "<div class=\"row text-center\">";
  foreach my $k (@{$tab})
  {
    if ($k =~ /(.*)-(.*)\.png$/)
    {
      $txt = $1.'-'.$2.'.txt';
      $name = $2;
    }
    print $h "</div><div class=\"row text-center\">" if $cmp != 0 && $cmp % 2 == 0;
    print $h "
    
    <div class=\"span6\">
    <h2>$name</h2>
    <p> <img src=\"$k\"/></p>
    <p class=\"featurette-p\"><a href=\"$txt\">text file</a></p>
    </div>
    ";
    $cmp++;
  }

  print $h "</div></div>";
}

sub mapnum
{
  my $dupmapnum = shift;
  my $dupnum_genome = shift;
  open (my $dupTE, $dupmapnum) || die "cannot open ".$dupmapnum."\n";
  my %dupnum_TE = ();
  my $header = <$dupTE>;
  while (<$dupTE>)
  {
    chomp $_;
    my @dupline = split /\t/, $_;
    $dupnum_TE{$dupline[0]} = $dupline[2];
  }
  close $dupTE;
  open (my $du_TE, '>'.$dupmapnum) || die "cannot open to write ".$dupmapnum."\n";
  print $du_TE "sequence\tduplicate\tgenome map num\tmap num\n";
  while (my ($k, $v) = each %dupnum_TE )
  {
    my $hashRef = ${$dupnum_genome}{$k};
    print $du_TE "$k\t$hashRef->[0]\t$hashRef->[1]\t$v\n";
  }
  close $du_TE;
}

sub copy_css
{
  my $dir = shift;
  my $path = dirname(__FILE__);
  dircopy( $path.'/css', $dir.'/css' );
}


sub copy_js
{
  my $dir = shift;
  my $path = dirname(__FILE__);
  dircopy( $path.'/js', $dir.'/js' );
}

1;

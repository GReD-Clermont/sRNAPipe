package sRNAPipe::Rcall;

use strict;
use warnings;
use Statistics::R;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( &histogram &pie_chart &bg_to_png );

sub histogram
{
    my ($size_hashR, $out_png, $size)    = @_;
    my (@abs, @ord);
    my $i = 0;
    foreach my $k (sort {$a <=> $b} keys %{$size_hashR})
    {
        my $percentage = 0;
        $percentage = $size_hashR->{$k} * 100 / $size if $size != 0;
        $abs[$i] = $k ; $ord[$i] = $percentage; $i++;
    }
    my $abs = join (",", @abs );
    my $ord = join (",", @ord );
    if (scalar(@abs) != 0)
    {

        my $R = Statistics::R->new();
        $R->startR;
        $R->send(
            qq`library(ggplot2)
            percentage = c($ord)
            size =c($abs)
            min = min(size)
            max = max(size)
            dat = data.frame(size,percentage)
            png(filename=\"$out_png\", width = 640, height = 640)
            c = ggplot(dat,aes(size,percentage))
            c + geom_bar(stat="identity") + scale_x_continuous(breaks=min:max)+theme( axis.text.x = element_text(angle=90, hjust=0.5, size=20), axis.text.y = element_text( size=20 ), axis.title.x = element_text( size=25, face="bold"), axis.title.y = element_text( size=25, face="bold") )
            dev.off()`);
        $R->stopR();

    }
}

sub bg_to_png
{
    my ( $fai, $bgP, $bgM, $dir, $sb ) = @_;
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(
    qq`library('Sushi')
    fai =read.table("$fai")
    if ( file.info("$bgP")\$size !=0 )
    {
        bgP = read.table("$bgP")
    } else { bgP = data.frame(factor(),integer()) }

    if ( file.info("$bgM")\$size !=0 )
    {
        bgM = read.table("$bgM")
    } else { bgM = data.frame(factor(),integer()) }

    f_both = function(chr,end) {
        jpeg( paste0("$dir",as.character(chr),".png"), quality=100)
        par(mfrow=c(2,1),mar=c(1,10,1,3))
        plotBedgraph(bgP, chrom=chr,chromstart=0,chromend=end,transparency=.50,    color=SushiColors(2)(2)[1])
        axis(side=2,las=2,tcl=.2)
        mtext("Scaled Read Depth",side=2,line=4,cex=1,font=2)
        plotBedgraph(bgM, chrom=chr,chromstart=0,chromend=end,transparency=.50, flip=TRUE, color=SushiColors(2)(2)[2])
        labelgenome(chrom=chr,chromstart=0,chromend=end,side=3,n=3,scale="$sb", line=0,    chromline = 0.5,    scaleline = 0.5, scaleadjust =1.05, chromadjust = -0.4)
        axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),labels=-1*pretty(par("yaxp")[c(1,2)]))
        mtext("Scaled Read Depth",side=2,line=4.5,cex=1,font=2)
        dev.off()
    }

    f_plus = function(chr,end) {
        jpeg( paste0("$dir",as.character(chr),".png"), quality=100)
        plotBedgraph(bgP, chrom=chr,chromstart=0,chromend=end,transparency=.50,    color=SushiColors(2)(2)[1])
        labelgenome(chrom=chr,chromstart=0,chromend=end,n=3,scale="$sb", line=0,    chromline = 0.5,    scaleline = 0.5, scaleadjust =1.05, chromadjust = -0.4)
        axis(side=2,las=2,tcl=.2)
        mtext("Scaled Read Depth",side=2,line=4,cex=1,font=2)
        dev.off()
    }

    f_minus = function(chr,end) {
        jpeg( paste0("$dir",as.character(chr),".png"), quality=100)
        plotBedgraph(bgM, chrom=chr,chromstart=0,chromend=end,transparency=.50, flip=TRUE, color=SushiColors(2)(2)[2])
        labelgenome(chrom=chr,chromstart=0,chromend=end,n=3,scale="$sb", line=0,    chromline = 0.5,    scaleline = 0.5, scaleadjust =1.05, chromadjust = -0.4)
        axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),labels=-1*pretty(par("yaxp")[c(1,2)]))
        mtext("Scaled Read Depth",side=2,line=4.5,cex=1,font=2)
        dev.off()
    }

    fai_b = fai[fai\$V1 %in% intersect(bgM\$V1,bgP\$V1), ]
    mapply( f_both, fai_b\$V1, fai_b\$V2)

    fai_p = fai[fai\$V1 %in% setdiff(bgP\$V1,bgM\$V1), ]
    mapply( f_plus, fai_p\$V1, fai_p\$V2)

    fai_m = fai[fai\$V1 %in% setdiff(bgM\$V1,bgP\$V1), ]
    mapply( f_minus, fai_m\$V1, fai_m\$V2) `);

    $R->stopR();
}

sub pie_chart
{
    my $dir = shift;
    my $in = $dir.'repartition.txt';
    my $out = $dir.'pie_chart.png';

    my $R = Statistics::R->new();
    $R->startR;
    $R->send(
    qq`
    library(plotrix)
    library(RColorBrewer)
    R =read.table("$in",header=T)
    values = round(R\$percentage)
    keys = R\$type
    lab = paste(values, "%", sep="")
    png("$out")
    colors <- brewer.pal(7,"Paired")
    pie(values, col=colors, labels=lab, clockwise=TRUE)
    legend("bottom", legend = keys, fill=colors, bty="n", ncol = 3)
    par(mai = c(0,0,0,0))
    layout(c(1,2),heights=c(0.3,1))
    plot.new()
    legend("bottom", legend = keys, fill=colors, bty="n",ncol = 3)
    pie(values, col=colors, labels=lab, clockwise=TRUE)
    dev.off()`
    );
    $R->stopR();
}

1;

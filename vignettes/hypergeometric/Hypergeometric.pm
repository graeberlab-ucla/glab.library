package Hypergeometric;

use strict;
use lib '/data6/jawong/MMHCC/perlscripts/modules';
#use lib '/home/software/perl5/site_perl/5.8.0/i386-linux-thread-multi'; #for Math:Libm on the cluster (and on rover and hydra)
use lib '/home/software/site_perl/5.8.0/i386-linux-thread-multi'; #for Math:Libm on the cluster (and on rover and hydra)
use Stats_tgg;
use Math::Libm; # on rover -> /usr/lib/perl5/site_perl/5.8.5/i386-linux-thread-multi/



use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(
	hypergeom
	binomialZ
	binomialZ_bidirectional
	kolmogorov_smirnov 
);


# get hypergeometric distribution; taken from hypergeometric.pl
sub hypergeom {

my ($nredtotal, $nredpick, $ntotal, $ntries) = @_;

my $result;
my $c3 = &binomial($ntotal,$ntries);

#print "$c3\n";

my $xmax=$ntries;
if ($nredtotal < $ntries) {$xmax = $nredtotal;}

for (my $x=$nredpick; $x<=$xmax; $x++) {

  my $c1 = &binomial($nredtotal,$x);
  my $c2 = &binomial($ntotal-$nredtotal,$ntries-$x);

  die "c1 $c1\n($nredtotal,$x)" if ($c1 =~ /ERROR/);
  die "c2 $c2\n($ntotal-$nredtotal,$ntries-$x)" if ($c2 =~ /ERROR/);

#    print "$c1\t$c2\n";

  my $y = $c1*$c2/$c3;
#  print "$y\t";
  $result += $y;
}

#print "\n$result\n";

return($result);

}


# approxmiates hypergeometric cumulative distribution using binomial approximated by using the normal distribution.
# taken from hypergeometric.pl and modified slightly
sub binomialZ { #in the limit N,M -> infinity
#see http://www.math.uah.edu/stat/bernoulli/index.html
#or /data1/users/graeber/reference/mathematics/stats/bernoulli_binomial_variance/index.html

    my ($nredtotal, $nredpick, $ntotal, $ntries) = @_;
    my $rtrue = $nredtotal/$ntotal;
    my $nexpect = $rtrue * $ntries;
    my $nsigma = sqrt($rtrue * (1-$rtrue) * $ntries);
    my $Z = ($nredpick - $nexpect)/$nsigma;
    #my $P = &erfc($Z);
	#my $P = 1 - 0.5 * &Math::Libm::erfc($Z/sqrt(2));		# cumulative distribution for Z: 0->z
	my $P = 0.5 * &Math::Libm::erfc($Z/sqrt(2));			# cumulative distribution for Z: z->inf
    return $P;
}


sub binomialZ_bidirectional { #in the limit N,M -> infinity
#see http://www.math.uah.edu/stat/bernoulli/index.html
#or /data1/users/graeber/reference/mathematics/stats/bernoulli_binomial_variance/index.html

    my ($nredtotal, $nredpick, $ntotal, $ntries) = @_;
    my $rtrue = $nredtotal/$ntotal;
    my $nexpect = $rtrue * $ntries;
    my $nsigma = sqrt($rtrue * (1-$rtrue) * $ntries);
    my $Z = abs($nredpick - $nexpect)/$nsigma;
    
    #my $P = &erfc($Z);
    my $P = 0.5 * &Math::Libm::erfc($Z/sqrt(2));
    #if using erfc from Math::Cephes or Math::Libm, need to divide the output by 2 and divide the Z input by sqrt 2, if using algorithms/erfc these steps are already incorporated
    
    $P = -1 * $P if ($nredpick - $nexpect < 0);
    
    if ($P eq "") {die "it appears the call to sub erfc is not working";}
    
    #print "Z $Z P $P\n";
    
    return $P;
}

sub kolmogorov_smirnov 
{

my ($nredtotal, $nredpick, $ntotal, $ntries) = @_;
#M m N n
#my ($samplesize, $samplesuccesses, $populationsize, $populationsuccesses) = @_;

# my $ks = $nredpick/$nredtotal - $ntries/$ntotal;
my $ks = $nredpick / $nredtotal - (($ntries - $nredpick) / ($ntotal - $nredtotal));
return($ks);

}

1;
__END__
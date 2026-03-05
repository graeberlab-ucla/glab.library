package Stats_tgg;
 
########################################################################
#
# Statistical analysis package
# 
#
#
#                                                     03/07/00  TGG
#
########################################################################
 
#add ~/perl/hypergeometric.pl


require Exporter;

@ISA    =qw(Exporter);

@EXPORT =qw(
            max
            min
            absmax
            absmin
            maxarray
            minarray
            absmaxarray
            absminarray
            mean
            meansd
            meansd0
            meansdblank
            meanvar
            meanvarblank
            meanvarskewblank
            meanvarexpblank
            median
            cofv
            cofvblank
            cofvsq
            cofvsqblank
            corrcoef
            corrcoefblank
            corrcoefpositiveblank
            corrcoef0
            corrcoef00
            corrcoef0blank
            corrcoefP
            corrcoefPblank
            corrcoefrndblank
            dotproduct
            dotproductunnorm
            dotproductblank
            dotproductPblank
            dotproductPPblank
            corrcoefPPblank
            corrcoefPPpositiveblank
            corrcoefP0
            corrcoefsubprof
            linear_regression0
            linear_regressionblank
            linear_regression_varianceblank
            linear_regression0blank
            principle0
            factorial
            factorial_ratio
            binomial
            binomial_sterling
            poisson
            significantdigits
            ttest
            ttest_w_error
            betai
            betacf
            gammln
            shufflearray
            arrayadd
            arraysub
            arrayave
            arrayaveblank
            matrixmult
            matrixadd
	    );
            #roc
            #avearray
# 'blank' versions eliminate component pairs when one from either vector is blank ("")
# '0' versions eliminate component pairs when both are zero
# '00' versions eliminate component pairs when one from either vector is zero
use strict;

#need one of the packages below for sub roc
#use Math::Libm ':all';  #erfc
#use Math::Cephes qw(:all);  #erfc
#use Math::Cephes qw(:misc);  #erfc

############################################################
#                                                          
# max
#                                                   
############################################################

sub max {
    my $max = shift;
    my $new = shift;
    if ($new > $max) {
        $max = $new;
    }
    return $max;
}


############################################################
#                                                          
# min
#                                                   
############################################################

sub min {
    my $min = shift;
    my $new = shift;
    if ($new < $min) {
        $min = $new;
    }
    return $min;
}


############################################################
#                                                          
# absmax
#                                                   
############################################################

sub absmax {
    my $max = shift;
    my $new = shift;
    if (abs($new) > abs($max)) {
        $max = $new;
    }
    return $max;
}


############################################################
#                                                          
# absmin
#                                                   
############################################################

sub absmin {
    my $min = shift;
    my $new = shift;
    if (abs($new) < abs($min)) {
        $min = $new;
    }
    return $min;
}


############################################################
#                                                          
# maxarray
#                                                   
############################################################

sub maxarray {
    my $point = shift;
    my $max = $point->[0];
    my @is = (0);
    for (my $i=1; $i<=$#{$point}; $i++) {
        if ($point->[$i] > $max) {
            $max = $point->[$i];
            @is = ($i);
        } elsif ($point->[$i] == $max) {
            push(@is, $i);
        }
    }
    return ($max, \@is);
}


############################################################
#                                                          
# minarray
#                                                   
############################################################

sub minarray {
    my $point = shift;
    my $min = $point->[0];
    my @is = (0);
    for (my $i=1; $i<=$#{$point}; $i++) {
        if ($point->[$i] < $min) {
            $min = $point->[$i];
            @is = ($i);
        } elsif ($point->[$i] == $min) {
            push(@is, $i);
        }
    }
    return ($min, \@is);
}


############################################################
#                                                          
# absmaxarray
#                                                   
############################################################

sub absmaxarray {
    my $point = shift;
    my $max = $point->[0];
    my @is = (0);
    for (my $i=1; $i<=$#{$point}; $i++) {
        if (abs($point->[$i]) > abs($max)) {
            $max = $point->[$i];
            @is = ($i);
        } elsif (abs($point->[$i]) == abs($max)) {
            push(@is, $i);
        }
    }
    return ($max, \@is);
}


############################################################
#                                                          
# absminarray
#                                                   
############################################################

sub absminarray {
    my $point = shift;
    my $min = $point->[0];
    my @is = (0);
    for (my $i=1; $i<=$#{$point}; $i++) {
        if (abs($point->[$i]) < abs($min)) {
            $min = $point->[$i];
            @is = ($i);
        } elsif (abs($point->[$i]) == abs($min)) {
            push(@is, $i);
        }
    }
    return ($min, \@is);
}


############################################################
#                                                          
# mean
#                                                   
############################################################

sub mean {

my $xpoint = shift;

my @x = @{$xpoint};

my ($error, $xi);
my ($xn, $xsum, $xbar);

foreach $xi (@x) {
    $xsum += $xi;
    $xn++;
}

if ($xn < 1) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;

return ($error, $xbar);

}


############################################################
#                                                          
# meansd - mean and standard deviation
#                                                   
############################################################

sub meansd {

my $xpoint = shift;

my @x = @{$xpoint};

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #die "$var" if ($var < -1e-13);
    die "$var" if ($var < -1e-11);
    if ($var !~ /\-3\.9790393202/) {  #a common value for one program
        #warn "variance of $var in sub meansd, assuming var = 0";
    }
    $sd=0;
} else {     
    $sd = sqrt($var);
}

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $x2sum, $xn);

}


############################################################
#                                                          
# meansd0 - mean and standard deviation ignoring zeros
#                                                   
############################################################

sub meansd0 {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] == 0);
     $x[$j] = $x_pre[$i];
     $j++;
}

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd);

}


############################################################
#                                                          
# meansdblank - mean and standard deviation ignoring blanks
#                                                   
############################################################

sub meansdblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $j++;
}

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 1) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return ($error, $xbar, "", $x2sum, $xn);
#    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return ($error, $xbar, "", $x2sum, $xn);
}

$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $x2sum, $xn);

}


############################################################
#                                                          
# meanvar - mean and standard deviation
#                                                   
############################################################

sub meanvar {

my $xpoint = shift;

my @x = @{$xpoint};

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #if ($var eq "-7.27595761418343e-12") {  #a common value for one program
    #    warn "variance of $var in sub meansd, assuming var = 0";
    #}
    #die "$var" if ($var < -1e-13);
    die "$var" if ($var < -1e-11);
    if ($var !~ /\-3\.9790393202/) {  #a common value for one program
        #warn "variance of $var in sub meanvar, assuming var = 0";
    }
    $var=0;
}

return ($error, $xbar, $var);

}


############################################################
#                                                          
# meanvarblank - mean and standard deviation ignoring blanks
#                                                   
############################################################

sub meanvarblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $j++;
}
# $j = number of elements non-blank

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #if ($var eq "-7.27595761418343e-12") {  #a common value for one program
    #    warn "variance of $var in sub meansd, assuming var = 0";
    #}
    #die "$var" if ($var < -1e-13);
    die "$var" if ($var < -1e-11);
    if ($var !~ /\-3\.9790393202/) {  #a common value for one program
        #warn "variance of $var in sub meanvar, assuming var = 0";
    }
    $var=0;
}

return ($error, $xbar, $var, $j); # $j = number of elements non-blank

}


############################################################
#                                                          
# meanvarskewblank - mean and standard deviation ignoring blanks
#                                                   
############################################################

sub meanvarskewblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $j++;
}
# $j = number of elements non-blank

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #if ($var eq "-7.27595761418343e-12") {  #a common value for one program
    #    warn "variance of $var in sub meansd, assuming var = 0";
    #}
    #die "$var" if ($var < -1e-13);
    die "$var" if ($var < -1e-11);
    if ($var !~ /\-3\.9790393202/) {  #a common value for one program
        #warn "variance of $var in sub meanvar, assuming var = 0";
    }
    $var=0;
}

my $skew;
#my ($ep, $adev, $var2);
foreach $xi (@x) {
    my $s = $xi - $xbar;
	#$ep = $ep + $s;
	#$adev = $adev + abs($s);
	my $p = $s * $s;
	#$var2 = $var2 + $p;
	$p = $p * $s; 
	$skew = $skew + $p; 
}

if ($var != 0) { 
	$skew = $skew/($xn*(sqrt($var)**3));
} else {
	$error = "undefined2";
    return $error;
}

return ($error, $xbar, $var, $skew, $j); # $j = number of elements non-blank

=comment
http://www.aip.de/groups/soe/local/numres/bookfpdf/f14-1.pdf 
numerical recipes
14.1 Moments of a Distribution: Mean, Variance, Skewness

s=0.
do j=1,n
	s=s+data(j) 
enddo 11

ave=s/n
adev=0.
var=0.
skew=0.
curt=0.
ep=0.

do
	j=1,n 
	s=data(j)-ave 
	ep=ep+s 
	adev=adev+abs(s) 
	p=s*s
	var=var+p
	p=p*s 
	skew=skew+p 
	p=p*s 
	curt=curt+p
enddo

adev=adev/n 
var=(var-ep**2/n)/(n-1) 
sdev=sqrt(var)

if(var.ne.0.)then 
	skew=skew/(n*sdev**3) 
	curt=curt/(n*var**2)-3. #kurtosis curtosis
=cut
}


############################################################
#                                                          
# meanvarexpblank - mean and standard deviation ignoring blanks
#                   using the exponential of all values
#                                                   
############################################################

sub meanvarexpblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = exp($x_pre[$i]);
     $j++;
}
# $j = number of elements non-blank

my ($error, $xi, $sd);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #if ($var eq "-7.27595761418343e-12") {  #a common value for one program
    #    warn "variance of $var in sub meansd, assuming var = 0";
    #}
    #die "$var" if ($var < -1e-13);
    die "$var" if ($var < -1e-11);
    if ($var !~ /\-3\.9790393202/) {  #a common value for one program
        #warn "variance of $var in sub meanvar, assuming var = 0";
    }
    $var=0;
}

return ($error, $xbar, $var, $j); # $j = number of elements non-blank

}


############################################################
#                                                          
# median
#
############################################################

sub median {
    my $x = shift;
    my @x = sort by_asc_median @{$x};
    
    my $n = $#x/2;
    my $nint = int($n);
    
    my $median;
    if ($n == $nint) {
        $median = $x[$n];
    } else {
        $median = ($x[$nint] + $x[$nint+1])/2;
    } 
    return $median;

    sub by_asc_median {
        $a <=> $b; 
    }
}



############################################################
#                                                          
# cofv - coefficient of variation, mean and standard deviation
#          cofv = sd/mean
#
############################################################

sub cofv {

my $xpoint = shift;

my @x = @{$xpoint};

my ($error, $xi, $sd, $cofv);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    die "$var" if ($var < -1e-14);
    warn "variance of $var in sub cofv, assuming var = 0";
    $sd=0;
} else {     
    #$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
    $sd = sqrt($var);
}

$cofv = "-";
if ($xbar != 0) {
    $cofv = $sd/$xbar;
}

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $cofv);

}


############################################################
#                                                          
# cofvblank - coefficient of variation, mean and standard deviation
#             ignoring blanks
#             cofv = sd/mean
#                                      
############################################################

sub cofvblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $j++;
}

my ($error, $xi, $sd, $cofv);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    #die "$var" if ($var < -1e-14);
    warn "variance of $var in sub cofvblank, assuming var = 0";
    $sd=0;
} else {     
    $sd = sqrt($var);
    #$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
}

$cofv = "-";
if ($xbar != 0) {
    $cofv = $sd/$xbar;
}

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $cofv);

}


############################################################
#
# cofvsq - coefficient of variation squared, mean and standard deviation
#             cofv = sd/mean
#
############################################################

sub cofvsq {

my $xpoint = shift;

my @x = @{$xpoint};

my ($error, $xi, $sd, $cofvsq);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    die "$var" if ($var < -1e-14);
    warn "variance of $var in sub cofvsq, assuming var = 0";
    $sd=0;
} else {     
    #$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
    $sd = sqrt($var);
}

$cofvsq = "-";
if ($xbar != 0) {
    $cofvsq = $var/($xbar**2);
}

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $cofvsq);

}


############################################################
#                                                          
# cofvsqblank - coefficient of variation squared, mean and standard deviation
#             ignoring blanks
#             cofv = sd/mean
#                                      
############################################################

sub cofvsqblank {

my $xpoint = shift;

my @x_pre = @{$xpoint};

my @x;

my $j = 0;
for (my $i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $j++;
}

my ($error, $xi, $sd, $cofvsq);
my ($xn, $xsum, $x2sum, $xbar, $xbar2);

foreach $xi (@x) {
    $xsum += $xi;
    $x2sum += ($xi)**2;
    $xn++;
}

if ($xn < 2) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$xbar2 = ($xbar)**2;

my $var = ($x2sum - $xn*$xbar2)/($xn-1);
if ($var < 0) {
    die "$var" if ($var < -1e-14);
    warn "variance of $var in sub cofvsqblank, assuming var = 0";
    $sd=0;
} else {     
    #$sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
    $sd = sqrt($var);
}

$cofvsq = "-";
if ($xbar != 0) {
    $cofvsq = $var/($xbar**2);
}

#my $x2b;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#my $sdv1 = sqrt($x2b/($xn-1));
#print "sd $sd\tsdv1 $sdv1\n";

return ($error, $xbar, $sd, $cofvsq);

}


############################################################
#                                                          
# corrcoef - calculate Pearson correlation coefficient
#                                                   
############################################################

sub corrcoef {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 


for ($i=0; $i<=$#x; $i++) {
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));

#$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );
$rbottom = ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) );

if ($rbottom < 0) {
    $error = "undefined_rbottom<0";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$rbottom = sqrt($rbottom);

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

# inefficient method:
#my ($x2b, $y2b, $xs2b, $ys2b);
#my $xbar = $xsum/$n;
#my $ybar = $ysum/$n;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#$xs2b = sqrt($x2b);
#foreach $yi (@y) {
#    $y2b += ($yi - $ybar)**2;
#}
#$ys2b = sqrt($y2b);
#if ($xs2b == 0 || $ys2b == 0) {
#    $error = "undefined";
##    $error2 = "one or more vectors with no variation";
#    return($error);
##    return($error, $null, $error2);
#}
#my $rtopineff;
#for ($i = 0; $i <= $#x; $i++) {
#    $rtopineff += ($x[$i]-$xbar)*($y[$i]-$ybar);
#}
#my $rineff = $rtopineff/($xs2b*$ys2b);
#print "r $r        rineff $rineff\n";
#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $r);

}


############################################################
#                                                          
# corrcoefblank - calculate Pearson correlation coefficient
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                                   
############################################################

sub corrcoefblank {

my $xpoint = shift;
my $ypoint = shift;
my $printifdie = shift;

my ($error, $i, $rtop, $rbottom, $r);

if ($xpoint eq "" || $ypoint eq "") {
	if ($xpoint eq "") {
		$error .= "\$xpoint eq \"\"; ";
		#$error .= "\$xpoint eq null string; ";
	}
	if ($ypoint eq "") {
		$error .= "\$ypoint eq \"\"; ";
		#$error .= "\$ypoint eq null string; ";
	}
	return ($error);
}

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 


for ($i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
#print "xx " if ($x[$i] == 0 && $y[$i] == 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));

$rbottom = ( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );
if ($rbottom < 0) {
    warn "\$rbottom value of $rbottom converted to zero (n=$n)";
    $rbottom = 0;
    #die "$printifdie\n";};
}
$rbottom = sqrt( $rbottom );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

# inefficient method:
#my ($x2b, $y2b, $xs2b, $ys2b);
#my $xbar = $xsum/$n;
#my $ybar = $ysum/$n;
#foreach $xi (@x) {
#    $x2b += ($xi - $xbar)**2;
#}
#$xs2b = sqrt($x2b);
#foreach $yi (@y) {
#    $y2b += ($yi - $ybar)**2;
#}
#$ys2b = sqrt($y2b);
#if ($xs2b == 0 || $ys2b == 0) {
#    $error = "undefined";
##    $error2 = "one or more vectors with no variation";
#    return($error);
##    return($error, $null, $error2);
#}
#my $rtopineff;
#for ($i = 0; $i <= $#x; $i++) {
#    $rtopineff += ($x[$i]-$xbar)*($y[$i]-$ybar);
#}
#my $rineff = $rtopineff/($xs2b*$ys2b);
#print "r $r        rineff $rineff\n";
#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $r, $n);

}


############################################################
#                                                          
# corrcoefpositiveblank - calculate Pearson correlation coefficient
#                 ignoring components in which  
#                 either entry (from either vector) is blank or negative
#                                                   
############################################################

sub corrcoefpositiveblank {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 


for ($i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    next if ($x[$i] < 0 || $y[$i] < 0);
#print "xx " if ($x[$i] == 0 && $y[$i] == 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

return ($error, $r, $n);

}


############################################################
#                                                          
# corrcoef0 - calculate Pearson correlation coefficient
#             ignoring components in which both entries 
#             (from each vector) are zero                                      
#
############################################################

sub corrcoef0 {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

$n = 0;
for ($i=0; $i <= $#x; $i++) {
    next if ($x[$i] == 0 && $y[$i] == 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

return ($error, $r);

}


############################################################
#                                                          
# corrcoef00 - calculate Pearson correlation coefficient
#              ignoring components in which either entry 
#              (from the two vectors) is zero                                      
#
############################################################

sub corrcoef00 {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

$n = 0;
for ($i=0; $i <= $#x; $i++) {
    next if ($x[$i] == 0 || $y[$i] == 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

return ($error, $r, $n);

}


############################################################
#                                                          
# corrcoef0blank - calculate Pearson correlation coefficient
#             ignoring components in which both entries 
#             (from each vector) are zero or either entry is blank
#
############################################################

sub corrcoef0blank {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 

if ($#x != $#y) {
    $error = "vectors of different length";
    
    print "$#x != $#y\n";
    
    return $error;
}

$n = 0;
for ($i=0; $i <= $#x; $i++) {
    next if ($x[$i] == 0 && $y[$i] == 0);
#    next if ($x[$i] > 2100);
    next if ($x[$i] eq "" || $y[$i] eq "");
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

return ($error, $r);

}


############################################################
#                                                          
# corrcoefP - calculate Pearson correlation coefficient
#             and P value by shuffling
#                                       
############################################################

sub corrcoefP {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;

if ($nrnd < 2) {
    my $error = "must use nrnd >= 2";
    return ($error);
}

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($i, $rtop, $rbottom, $r);

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 

for ($i=0; $i<=$#x; $i++) {
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$xsumysum = $xsum*$ysum;
$rtop = (($n*$xysum)-($xsumysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
#    $error2 = "one or more vectors with no variation";
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;


#calculate P value
my @rrs;
for ($i=0; $i<$nrnd; $i++) {

    my $xysumr;

    #generate random vector
    my @ind;
    for (my $j=0; $j<=$n-1; $j++) {
        $ind[$j]=$j;
    }

    for (my $j=0; $j<=$n-1; $j++) {
        my $l  = $n - $j;
        my $rnd = int(rand $l);
        my $newj = splice(@ind, $rnd, 1);
        my $xrcomp = $x[$newj];
        $xysumr += $xrcomp*$y[$j];
    }
    #die if ($#ind != -1); #debugging

    my $rtopr = (($n*$xysumr)-($xsumysum));
    my $rr = $rtopr/$rbottom;

    push (@rrs, $rr);
}

my ($rrmean, $rrsd);
{
    my @x = @rrs;
    
    my ($xi, $sd);
    my ($xn, $xsum, $x2sum, $xbar, $xbar2);
    
    foreach $xi (@x) {
        $xsum += $xi;
        $x2sum += ($xi)**2;
        $xn++;
    }
    
    if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP: undefined:  nrnd ($nrnd) must be greater than 2";
        return $error;
    }
    
    $xbar = $xsum/$xn;
    $xbar2 = ($xbar)**2;
    
    $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
    
    ($rrmean, $rrsd) = ($xbar, $sd);
}

my $Z;
my $error;
if ($rrsd != 0) {
    $Z = abs(($r-$rrmean)/$rrsd);
} else {
    $error = "Z undefined";
    $Z = -999;
}

return ($error, $r, $Z, \@rrs);

}


############################################################
#                                                          
# corrcoefPblank - calculate Pearson correlation coefficient
#             and Z score by shuffling
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                       
############################################################

sub corrcoefPblank {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
my ($xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


#calculate Z score
my $Z;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    
    for (my $i=0; $i<$nrnd; $i++) {
    
        my ($nn, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq ""); #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            $xsum += $xr;
            $ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        die "why nn = 0 ?" if ($nn == 0);
        
        $rtop = (($nn*$xysum)-($xsum*$ysum));
        my $xvar = ($nn*$x2sum)-($xsum**2);
        my $yvar = ($nn*$y2sum)-($ysum**2);
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $rr;
        
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
        }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefP0 subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
}
return ($error, $r, $Z);

}
############################################################
#                                                          
# corrcoefrndblank - calculate Pearson correlation coefficient
#             between array 1 and shuffled array 2
#                                       
############################################################

sub corrcoefrndblank {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my ($norig);

$norig = $#x + 1;

if ($norig == 0) {
    my $error = "undefined";
    return($error);
}

    my $var0count=0;
    my $rr;
    
    for (my $i=0; $i<=0; $i++) {
    
        my ($nn, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq ""); #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            $xsum += $xr;
            $ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        die "why nn = 0 ?" if ($nn == 0);
        
        $rtop = (($nn*$xysum)-($xsum*$ysum));
        my $xvar = ($nn*$x2sum)-($xsum**2);
        my $yvar = ($nn*$y2sum)-($ysum**2);
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } 
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count in corrcoefP0 subroutine";
    }
   
my $error;    
return ($error, $rr);

}


############################################################
#                                                          
# dotproduct - calculate dotproduct
#  
# this actually calculates cos(angle between vectors)
# use routine dotproductunnorm for standard dotproduct
# need to fix names eventually
#                                     
############################################################

sub dotproduct {

my $xpoint = shift;
my $ypoint = shift;

my $error; 

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
#my ($xsum, $ysum);
my ($x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    #$xsum += $x[$i];
    #$ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = $xysum;
$rbottom = sqrt( ( $x2sum )*( $y2sum ) );

if ($rbottom == 0) {
    $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


return ($error, $r);

}


############################################################
#                                                          
# dotproductunnorm - calculate dotproduct
# 
# dotproduct subroutine above is actually cos(angle between vectors)
# need to fix names eventually
#                                       
############################################################

sub dotproductunnorm {

my $xpoint = shift;
my $ypoint = shift;

my $error; 

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my ($n);
my ($xysum); 

for (my $i=0; $i<=$#x; $i++) {
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

return ($error, $xysum);

}


############################################################
#                                                          
# dotproductblank - calculate dotproduct
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                       
# this actually calculates cos(angle between vectors)
# use routine dotproductunnorm* for standard dotproduct
# need to fix names eventually
#                                     
############################################################

sub dotproductblank {

my $xpoint = shift;
my $ypoint = shift;

my $error; 

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
#my ($xsum, $ysum);
my ($x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    #$xsum += $x[$i];
    #$ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = $xysum;
$rbottom = sqrt( ( $x2sum )*( $y2sum ) );

if ($rbottom == 0) {
    $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


return ($error, $r);

}


############################################################
#                                                          
# dotproductPblank - calculate dot product
#             and Z score by shuffling
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                       
# this actually calculates cos(angle between vectors)
# use routine dotproductunnorm* for standard dotproduct
# need to fix names eventually
#                                     
############################################################

sub dotproductPblank {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
#my ($xsum, $ysum);
my ($x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    #$xsum += $x[$i];
    #$ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = $xysum;
$rbottom = sqrt( ( $x2sum )*( $y2sum ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


#calculate Z score
my $Z;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    
    for (my $i=0; $i<$nrnd; $i++) {
        
        #my ($xsum, $ysum);
        my ($nn, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq ""); #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            #$xsum += $xr;
            #$ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        die "why nn = 0 ?" if ($nn == 0);
        
        $rtop = $xysum;
        my $xvar = $x2sum;
        my $yvar = $y2sum;
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $rr;
        
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
        }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefP0 subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        #my $sdsquared = (($x2sum - $xn*$xbar2)/($xn-1));
        #if ($sdsquared < 0) {
        #    my $old = $sdsquared;
        #    $sdsquared = abs($sdsquared);
        #    warn "converted $old to positive $sdsquared";
        #}

        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        #$sd = sqrt($sdsquared);
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
}
return ($error, $r, $Z);

}


############################################################
#                                                          
# dotproductPPblank - calculate dot product
#             and Z score, and P value by shuffling
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                       
# this actually calculates cos(angle between vectors)
# use routine dotproductunnorm* for standard dotproduct
# need to fix names eventually
#                                     
############################################################

sub dotproductPPblank {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
#my ($xsum, $ysum);
my ($x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    #$xsum += $x[$i];
    #$ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = $xysum;
$rbottom = sqrt( ( $x2sum )*( $y2sum ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


#calculate Z score, and P value
my $Z;
my $P;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    my $ngt;
    
    for (my $i=0; $i<$nrnd; $i++) {
        
        #my ($xsum, $ysum);
        my ($nn, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq ""); #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            #$xsum += $xr;
            #$ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        die "why nn = 0 ?" if ($nn == 0);
        
        $rtop = $xysum;
        my $xvar = $x2sum;
        my $yvar = $y2sum;
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $rr;
        
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
            if (abs($rr) > abs($r)) {
                $ngt++;
            }
        }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefP0 subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        #my $sdsquared = (($x2sum - $xn*$xbar2)/($xn-1));
        #if ($sdsquared < 0) {
        #    my $old = $sdsquared;
        #    $sdsquared = abs($sdsquared);
        #    warn "converted $old to positive $sdsquared";
        #}

        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        #$sd = sqrt($sdsquared);
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
    $P = $ngt/$nrnd;

}
return ($error, $r, $Z, $P);

}


############################################################
#                                                          
# corrcoefPPblank - calculate Pearson correlation coefficient
#             Z score, and P value by shuffling
#                 ignoring components in which  
#                 either entry (from either vector) is blank
#                                       
############################################################

sub corrcoefPPblank {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
my ($xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


#calculate Z score, and P value
my $Z;
my $P;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    my $ngt;
    
    for (my $i=0; $i<$nrnd; $i++) {
    
        my ($nn, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq ""); #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            $xsum += $xr;
            $ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        warn "why random nn = 0 ?" if ($nn == 0);
        
        $rtop = (($nn*$xysum)-($xsum*$ysum));
        my $xvar = ($nn*$x2sum)-($xsum**2);
        my $yvar = ($nn*$y2sum)-($ysum**2);
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $rr;
        
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
            if (abs($rr) > abs($r)) {
                $ngt++;
            }
        }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefP0 subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
    $P = $ngt/$nrnd;

}
return ($error, $r, $Z, $P);

}


############################################################
#                                                          
# corrcoefPPpositiveblank - calculate Pearson correlation coefficient
#             Z score, and P value by shuffling
#                 ignoring components in which  
#                 either entry (from either vector) is blank or negative
#                                       
############################################################

sub corrcoefPPpositiveblank {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
my ($xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");
    next if ($x[$i] < 0 || $y[$i] < 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;

my $nmin4Z = 4; 
if ($n <= $nmin4Z) {
    my $error = "\$n ($n) less than \$nmin4Z ($nmin4Z)";
    my $Z = "";
    my $P = "";
    return ($error, $r, $Z, $P, $n);
}

#calculate Z score, and P value
my $Z;
my $P;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    my $ngt;
    
    for (my $i=0; $i<$nrnd; $i++) {
    
        my ($nn, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            #next if ($xr eq "" && $y eq "");  #changed to below on 1/25/01
            next if ($xr eq "" || $y eq "");
            next if ($xr < 0 || $y < 0);
            $xsum += $xr;
            $ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        #######die "why nn = 0 ?" if ($nn == 0);
        
        my $rr;
       if ($nn==0) {
            $rr=0;
       } else {
        $rtop = (($nn*$xysum)-($xsum*$ysum));
        my $xvar = ($nn*$x2sum)-($xsum**2);
        my $yvar = ($nn*$y2sum)-($ysum**2);
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
            if (abs($rr) > abs($r)) {
                $ngt++;
            }
        }
       }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefPPpositiveblank subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
    $P = $ngt/$nrnd;

}
return ($error, $r, $Z, $P, $n);

}


############################################################
#                                                          
# corrcoefP0 - calculate Pearson correlation coefficient
#             and P value by shuffling
#             ignoring components in which both entries 
#             (from each vector) are zero
#                                       
############################################################

sub corrcoefP0 {

my $xpoint = shift;
my $ypoint = shift;
my $nrnd = shift;
my $ccmin = shift;
my $anticcmin = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

if ($#x != $#y) {
    my $error = "vectors of different length";
    return $error;
}

my $r;
my ($n, $norig);
my ($xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
my ($rtop, $rbottom);

$norig = $#x + 1;

for (my $i=0; $i<=$#x; $i++) {
    next if ($x[$i] == 0 && $y[$i] == 0);
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    my $error = "undefined";
    return($error);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    my $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;


#calculate Z score
my $Z;
my $error;
if ($r > $ccmin || $r < -$anticcmin) {
    #print "started a Z...";
    
    my @rrs;
    my $var0count=0;
    
    for (my $i=0; $i<$nrnd; $i++) {
    
        my ($nn, $xsum, $ysum, $x2sum, $y2sum, $xysum, $xsumysum); 
        my ($rtop, $rbottom);
        my @ind;
        for (my $j=0; $j<=$norig-1; $j++) {
            $ind[$j]=$j;
        }
    
        for (my $j=0; $j<=$norig-1; $j++) {
            my $l  = $norig - $j;
            my $rnd = int(rand $l);
            my $newj = splice(@ind, $rnd, 1);
            my $xr = $x[$newj];
            my $y = $y[$j];
            next if ($xr == 0 && $y == 0);
            $xsum += $xr;
            $ysum += $y;
            $x2sum += ($xr)**2;
            $y2sum += ($y)**2;
            $xysum += $xr*$y;
            $nn++;
        }
        #die if ($#ind != -1);  #debugging
    
        die "why nn = 0 ?" if ($nn == 0);
        
        $rtop = (($nn*$xysum)-($xsum*$ysum));
        my $xvar = ($nn*$x2sum)-($xsum**2);
        my $yvar = ($nn*$y2sum)-($ysum**2);
        $rbottom = sqrt( ( $xvar )*( $yvar ) );
    
        my $rr;
        
        my $var0;
        if ($rbottom == 0) {
            if ($xvar == 0 && $yvar == 0) {
                $rr = 1;
            } else {
                $var0=1;
                #die "why xbottom = 0 ?";
            }
        } else {
            $rr = $rtop/$rbottom;
        }
    
        if ($var0 == 1) {
            $var0count++;
            $i--;
        } else {
            push (@rrs, $rr);
        }
    }
    if ($var0count != 0) {
        warn "WARNING:  var0count $var0count (nrnd $nrnd) in corrcoefP0 subroutine";
    }
    
    #foreach my $rr (@rrs) {
    #    print $rr,"\n";
    #}
        
    my ($rrmean, $rrsd);
    {
        my @x = @rrs;
        
        my ($xi, $sd);
        my ($xn, $xsum, $x2sum, $xbar, $xbar2);
        
        foreach $xi (@x) {
            $xsum += $xi;
            $x2sum += ($xi)**2;
            $xn++;
        }
        
        if ($xn < 2) {
        my $error = "Stats_tgg.pm::corrcoefP0: undefined:  nrnd ($nrnd) must be greater than 2";
            return $error;
        }
        
        $xbar = $xsum/$xn;
        $xbar2 = ($xbar)**2;
        
        $sd = sqrt(($x2sum - $xn*$xbar2)/($xn-1));
        
        ($rrmean, $rrsd) = ($xbar, $sd);
    }    
    
    if ($rrsd != 0) {
        $Z = abs(($r-$rrmean)/$rrsd);
    } else {
        $error = "Z undefined";
        $Z = -999;
    }
#} else { # !($r > $ccmin || $r < -$anticcmin)
#    $Z = "";
}
return ($error, $r, $Z);

}


#sub corrcoefP00 {
#}


############################################################
#                                                          
# corrcoefsubprof - calculate Pearson correlation coefficient
#                                                   
############################################################

sub corrcoefsubprof {

my $xpoint = shift;
my $ypoint = shift;
my $dellibs = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $rtop, $rbottom, $r);

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

my %dellibsh;
foreach my $lib (@{$dellibs}) {
    $dellibsh{$lib}=1;
}

my ($n, $xsum, $ysum, $x2sum, $y2sum, $xysum); 


for ($i=0; $i<=$#x; $i++) {
    next if (exists $dellibsh{$i});
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

if ($n == 0) {
    $error = "undefined";
    return($error);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
    return($error);
}

$r = $rtop/$rbottom;

return ($error, $r);

}


############################################################
#                                                          
# linear_regression0 - calculate linear regression fit
#             ignoring components in which both entries 
#             (from each vector) are zero                                      
#
#could make this code more efficient (see corrcoef or linear_regressionblank vs linear_regressionblank_lessefficient)
#
############################################################

sub linear_regression0 {

my $xpoint = shift;
my $ypoint = shift;

my @x_pre = @{$xpoint};
my @y_pre = @{$ypoint};

my (@x, @y);

my ($error, $i, $j, $rtop, $r, $xi, $yi);

my ($xn, $yn, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b, $xysum, $x2sum);


if ($#x_pre != $#y_pre) {
    $error = "vectors of different length";
    return $error;
}

$j = 0;
for ($i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] == 0 && $y_pre[$i] == 0);
     $x[$j] = $x_pre[$i];
     $y[$j] = $y_pre[$i];
     $j++;
}

foreach $xi (@x) {
    $xsum += $xi;
    $xn++;
}

foreach $yi (@y) {
    $ysum += $yi;
    $yn++;
}

die if ($xn != $yn || $#x != $#y);
my $n = $xn;

if ($n == 0) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$n;
$ybar = $ysum/$n;

foreach $xi (@x) {
    $x2b += ($xi - $xbar)**2;
}
$xs2b = sqrt($x2b);

foreach $yi (@y) {
    $y2b += ($yi - $ybar)**2;
}
$ys2b = sqrt($y2b);

if ($xs2b == 0 || $ys2b == 0) {
    $error = "undefined";
#    $error = "one or more vectors with no variation";
    return $error;
}

for ($i = 0; $i <= $#x; $i++) {
    $rtop += ($x[$i]-$xbar)*($y[$i]-$ybar);
    $xysum += $x[$i]*$y[$i];
    $x2sum += $x[$i]**2;
}

$r = $rtop/($xs2b*$ys2b);

# y = ax + b
my $a = (($n*$xysum) - ($xsum*$ysum))/(($n*$x2sum) - $xsum**2);
my $b = $ybar - ($a*$xbar);

#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $a, $b, $r, $xbar, $ybar);

}


############################################################
#                                                          
# linear_regressionblank - calculate linear regression fit
#             ignoring components in which either entry is blank
#
#
############################################################

sub linear_regressionblank {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $j, $rtop, $rbottom, $r, $xi, $yi);

my ($n, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b, $xysum, $x2sum, $y2sum);

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

$n = 0;
for ($i=0; $i <= $#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");  #don't use components where either vector's value is blank
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

$xbar = $xsum/$n;
$ybar = $ysum/$n;

# y = ax + b
my $a = (($n*$xysum) - ($xsum*$ysum))/(($n*$x2sum) - $xsum**2);
my $b = $ybar - ($a*$xbar);

#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $a, $b, $r, $xbar, $ybar);

}


############################################################
#                                                          
# linear_regression_varianceblank - calculate linear regression fit
#             ignoring components in which either entry is blank
#
#
############################################################

sub linear_regression_varianceblank {

my $xpoint = shift;
my $ypoint = shift;

my @x = @{$xpoint};
my @y = @{$ypoint};

my ($error, $i, $j, $rtop, $rbottom, $r, $xi, $yi);

my ($n, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b, $xysum, $x2sum, $y2sum);

if ($#x != $#y) {
    $error = "vectors of different length";
    return $error;
}

$n = 0;
for ($i=0; $i <= $#x; $i++) {
    next if ($x[$i] eq "" || $y[$i] eq "");  #don't use components where either vector's value is blank
    $xsum += $x[$i];
    $ysum += $y[$i];
    $x2sum += ($x[$i])**2;
    $y2sum += ($y[$i])**2;
    $xysum += $x[$i]*$y[$i];
    $n++;
}

#my $error2;  # 9.20.00 $error2 added for debugging, not used by programs in general
#my $null;
if ($n == 0) {
    $error = "undefined";
#    $error2 = "one or more null vectors";
    return($error);
#    return($error, $null, $error2);
}

$rtop = (($n*$xysum)-($xsum*$ysum));
$rbottom = sqrt( ( ($n*$x2sum)-($xsum**2) )*( ($n*$y2sum)-($ysum**2) ) );

if ($rbottom == 0) {
    $error = "undefined";
#    $error2 = "one or more vectors with no variation";
    return($error);
#    return($error, $null, $error2);
}

$r = $rtop/$rbottom;

$xbar = $xsum/$n;
$ybar = $ysum/$n;

# y = ax + b
my $a = (($n*$xysum) - ($xsum*$ysum))/(($n*$x2sum) - $xsum**2);
my $b = $ybar - ($a*$xbar);

if (0) {
 my $b2 = ($x2sum*$ysum - $xsum*$xysum)/(($n*$x2sum) - $xsum**2);
 if ($b != $b2) {print "diff?\t$b\t$b2\n";}
}
#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

my $yvar = "";  #calculate residual
if ($n <= 2) {
    $error = "two few measurements to calculated variance in y => undetermined";
    return ($error, $a, $b, $r, $xbar, $ybar, $yvar);
} else {
    #reference: Taylor, J.R., An Introduction to Error Analysis, p. 158
    for ($i=0; $i <= $#x; $i++) {
        next if ($x[$i] eq "" || $y[$i] eq "");  #don't use components where either vector's value is blank
        $yvar += ($y[$i] - $a*$x[$i] - $b)**2;
    }
    $yvar = $yvar/($n-2);
}
return ($error, $a, $b, $r, $xbar, $ybar, $yvar);

}


############################################################
#                                                          
# linear_regressionblank_lessefficient - calculate linear regression fit
#             ignoring components in which either entry is blank
#
# use linear_regressionblank instead
#
############################################################

sub linear_regressionblank_lessefficient {
die "use linear_regressionblank instead";

my $xpoint = shift;
my $ypoint = shift;

my @x_pre = @{$xpoint};
my @y_pre = @{$ypoint};

my (@x, @y);

my ($error, $i, $j, $rtop, $r, $xi, $yi);

my ($xn, $yn, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b, $xysum, $x2sum);


if ($#x_pre != $#y_pre) {
    $error = "vectors of different length";
    return $error;
}

$j = 0;
for ($i=0; $i <= $#x_pre; $i++) {
     #next if ($x_pre[$i] == 0 && $y_pre[$i] == 0);
     next if ($x_pre[$i] eq "" || $y_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $y[$j] = $y_pre[$i];
     $j++;
}

foreach $xi (@x) {
    $xsum += $xi;
    $xn++;
}

foreach $yi (@y) {
    $ysum += $yi;
    $yn++;
}

die if ($xn != $yn || $#x != $#y);
my $n = $xn;

if ($n == 0) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$n;
$ybar = $ysum/$n;

foreach $xi (@x) {
    $x2b += ($xi - $xbar)**2;
}
$xs2b = sqrt($x2b);

foreach $yi (@y) {
    $y2b += ($yi - $ybar)**2;
}
$ys2b = sqrt($y2b);

if ($xs2b == 0 || $ys2b == 0) {
    $error = "undefined";
#    $error = "one or more vectors with no variation";
    return $error;
}

for ($i = 0; $i <= $#x; $i++) {
    $rtop += ($x[$i]-$xbar)*($y[$i]-$ybar);
    $xysum += $x[$i]*$y[$i];
    $x2sum += $x[$i]**2;
}

$r = $rtop/($xs2b*$ys2b);

# y = ax + b
my $a = (($n*$xysum) - ($xsum*$ysum))/(($n*$x2sum) - $xsum**2);
my $b = $ybar - ($a*$xbar);

#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $a, $b, $r, $xbar, $ybar);

}


############################################################
#                                                          
# linear_regression0blank - calculate linear regression fit
#             ignoring components in which both entries 
#             (from each vector) are zero or either entry is blank
#
#could make this code more efficient (see corrcoef or linear_regressionblank vs linear_regressionblank_lessefficient)
#
############################################################

sub linear_regression0blank {

my $xpoint = shift;
my $ypoint = shift;

my @x_pre = @{$xpoint};
my @y_pre = @{$ypoint};

my (@x, @y);

my ($error, $i, $j, $rtop, $r, $xi, $yi);

my ($xn, $yn, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b, $xysum, $x2sum);


if ($#x_pre != $#y_pre) {
    $error = "vectors of different length";
    return $error;
}

$j = 0;
for ($i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] == 0 && $y_pre[$i] == 0);
     next if ($x_pre[$i] eq "" || $y_pre[$i] eq "");
     $x[$j] = $x_pre[$i];
     $y[$j] = $y_pre[$i];
     $j++;
}

foreach $xi (@x) {
    $xsum += $xi;
    $xn++;
}

foreach $yi (@y) {
    $ysum += $yi;
    $yn++;
}

die if ($xn != $yn || $#x != $#y);
my $n = $xn;

if ($n == 0) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$n;
$ybar = $ysum/$n;

foreach $xi (@x) {
    $x2b += ($xi - $xbar)**2;
}
$xs2b = sqrt($x2b);

foreach $yi (@y) {
    $y2b += ($yi - $ybar)**2;
}
$ys2b = sqrt($y2b);

if ($xs2b == 0 || $ys2b == 0) {
    $error = "undefined";
#    $error = "one or more vectors with no variation";
    return $error;
}

for ($i = 0; $i <= $#x; $i++) {
    $rtop += ($x[$i]-$xbar)*($y[$i]-$ybar);
    $xysum += $x[$i]*$y[$i];
    $x2sum += $x[$i]**2;
}

$r = $rtop/($xs2b*$ys2b);

# y = ax + b
my $a = (($n*$xysum) - ($xsum*$ysum))/(($n*$x2sum) - $xsum**2);
my $b = $ybar - ($a*$xbar);

#print "$xn\t$xsum\t$xbar\t$x2b\t$xs2b\n";
#print "$yn\t$ysum\t$ybar\t$y2b\t$ys2b\n";
#print "$rtop\t$r\n";
#print "$r\n";

return ($error, $a, $b, $r, $xbar, $ybar);

}


############################################################
#                                                          
# principle0 - calculate individual 'components' of Pearson 
#              correlation coefficient
#              ignoring vector components in which both entries 
#              (from each vector) are zero                                      
#
#could make this code more efficient (see corrcoef or linear_regressionblank vs linear_regressionblank_lessefficient)
#
############################################################

sub principle0 {

my $xpoint = shift;
my $ypoint = shift;

my @x_pre = @{$xpoint};
my @y_pre = @{$ypoint};

my (@x, @y, @rtop, @r);

my ($error, $i, $j, $rtop, $rbottom, $r, $xi, $yi);

my ($xn, $yn, $xsum, $ysum, $xbar, $ybar, $x2b, $y2b, $xs2b, $ys2b);

#@r=@rtop=();

if ($#x_pre != $#y_pre) {
    $error = "vectors of different length";
    return $error;
}

$j = 0;
for ($i=0; $i <= $#x_pre; $i++) {
     next if ($x_pre[$i] == 0 && $y_pre[$i] == 0);
     $x[$j] = $x_pre[$i];
     $y[$j] = $y_pre[$i];
     $j++;
}

foreach $xi (@x) {
    $xsum += $xi;
    $xn++;
}

foreach $yi (@y) {
    $ysum += $yi;
    $yn++;
}

if ($xn == 0 || $yn == 0) {
    $error = "undefined";
#    $error = "one or more null vectors";
    return $error;
}

$xbar = $xsum/$xn;
$ybar = $ysum/$yn;

foreach $xi (@x) {
    $x2b += ($xi - $xbar)**2;
}
$xs2b = sqrt($x2b);

foreach $yi (@y) {
    $y2b += ($yi - $ybar)**2;
}
$ys2b = sqrt($y2b);

if ($xs2b == 0 || $ys2b == 0) {
    $error = "undefined";
#    $error = "one or more vectors with no variation";
    return $error;
}

$rbottom = ($xs2b*$ys2b);
die if ($rbottom == 0);

for ($i = 0; $i <= $#x; $i++) {
    $rtop[$i] = ($x[$i]-$xbar)*($y[$i]-$ybar);
#    if ($rtop[$i] eq "") {$rtop[$i] = 0};    #doesn't do any good
    $r[$i] = $rtop[$i]/$rbottom;
#    if ($r[$i] eq "") {$r[$i] = 0};    
}

return ($error, \@rtop, \@r);

}


############################################################
#                                                          
# factorial        $k!
#
############################################################

sub factorial {
    my $k = shift;
    
    if ($k == 0) {return(1)};
    
    if (int $k != $k || $k < 0) {
        return("ERROR:input to factorial must be an integer greater than or equal to zero");
    }
    
    my $kf = 1;
    for (my $i=2;$i<=$k;$i++) {
        $kf *= $i;
    }
    
    return($kf);
     
}


############################################################
#                                                          
# factorial_ratio        $k1!/$k2!
#
############################################################

sub factorial_ratio {
    my $k1 = shift;
    my $k2 = shift;

    if (int $k1 != $k1 || $k1 < 0 || int $k2 != $k2 || $k2 < 0) {
        return("ERROR:inputs to factorial_ratio must be integers greater than or equal to zero");
    }
    
    if ($k1 - $k2 == 0) {return(1)};
    
    my $invert;
    if ($k1 < $k2) { # force $k1 to be greater than $k2
        $invert = 1;
        my $x = $k1;
        $k1 = $k2;
        $k2 = $x;
    } 
    
    if ($k2 == 0) {$k2 = 1};
        
    my $kfr = 1;
    for (my $i=$k2+1;$i<=$k1;$i++) {
        $kfr *= $i;
    }
        
    if ($invert) {
        $kfr = 1/$kfr;
    }
    
    return($kfr);
     
}


############################################################
#                                                          
# binomial        $k1! / ( ($k1-$k2)! * $k2! )
#
#     note:  binomial (n,m) = binomial (n,n-m)
#
############################################################

sub binomial {
    my $k1 = shift;
    my $k2 = shift;

    if (int $k1 != $k1 || $k1 < 0 || int $k2 != $k2 || $k2 < 0) {
        return("ERROR:inputs to binomial must be integers greater than or equal to zero");
    }
    
    if ($k1 < $k2) { # force $k1 to be greater than $k2
        return("ERROR:first input to binomial must be larger than second input");
    } 

    if ($k1 - $k2 == 0 || $k2 == 0) {return(1)};
    
    #determine whether $k1-$k2 or $k2 is larger
    my ($x, $y);
    if (($k1-$k2) > $k2) {
        $x = $k1-$k2;
        $y = $k2;
    } else {
        $x = $k2;
        $y = $k1-$k2;
    } #now x is larger than y, and $x + $y = $k1

    my $kfr = 1;
    for (my $i=1;$i<=$y;$i++) {
        $kfr *= ($i+$x)/$i;
    }

    return($kfr);

}


############################################################
#                                                          
# binomial_sterling        $k1! / ( ($k1-$k2)! * $k2! )
#
#     note:  binomial (n,m) = binomial (n,n-m)
#
#     using Sterling's apporimation  n! =~ (e**-n)(n**n)sqrt(2 pi n)
#
############################################################

sub binomial_sterling {
    my $k1 = shift;
    my $k2 = shift;

    my $sqrt2pi = sqrt(2 * 3.141592654);

    if (int $k1 != $k1 || $k1 < 0 || int $k2 != $k2 || $k2 < 0) {
        return("ERROR:inputs to binomial must be integers greater than or equal to zero");
    }
    
    if ($k1 < $k2) { # force $k1 to be greater than $k2
        return("ERROR:first input to binomial must be larger than second input");
    } 

    if ($k1 - $k2 == 0 || $k2 == 0) {return(1)};
    
    #determine whether $k1-$k2 or $k2 is larger
    my ($x, $y);
    if (($k1-$k2) > $k2) {
        $x = $k1-$k2;
        $y = $k2;
    } else {
        $x = $k2;
        $y = $k1-$k2;
    } #now x is larger than y, and $x + $y = $k1

    my $z = $x-$y;
    my $kfr = ( ($x**$x) / ( $sqrt2pi * ($y**$y) * ($z**$z) ) ) * (sqrt ($x/($y*$z)));

    return($kfr);

}


############################################################
#                                                          
# poisson        ( e**(-lamda) * (lamda**n) ) / n!
#
#
############################################################

sub poisson {
    my $m = shift; #lamda, mean
    my $n = shift;

    if (int $m != $m || $m < 0 || int $n != $n || $n < 0) {
        return("ERROR:inputs to poisson must be integers greater than or equal to zero");
    }
     
    my $numcount = $n;
    my $denomcount = $n; 

    my $x = 2.718281828**(-$m);
    
    while ($numcount || $denomcount) {
        if ($x >= 1 && $numcount) {
            $x *= $m;
            last if ($x =~ /INF|NaNQ/);
            $numcount--;
        } elsif ($x < 1 && $denomcount) {
            $x /= $denomcount;
            last if ($x =~ /INF|NaNQ/);
            $denomcount--;
        } elsif ($numcount) {
            $x *= $m;
            last if ($x =~ /INF|NaNQ/);
            $numcount--;
        } elsif ($denomcount) {
            $x /= $denomcount;
            last if ($x =~ /INF|NaNQ/);
            $denomcount--;
        }
    }
    
    return($x);

}


############################################################
#                                                          
# significantdigits
#
#
############################################################

sub significantdigits {
    my $x = shift;
    my $n = shift;
    
    die "sub significantdigits only coded for values between 0 and 1" if ($x > 1 || $x <0); 
    die "must want at least one significant digit" if ($n < 1);
    die "number of desired sig. digits must be an integer" if ($n != int($n));

    if ($x == 0) {return("0");}
    if ($x =~ /^1/) {
      if ($x =~ /^1\.{0,1}$/) {return("1");}
      if ($x =~ /^1\.(0+)$/) {
        my $nzeros = length($1);
        if ($n == 1) {return("1");}
        my $string = "1.";
        for (my $i=1; $i<=$n-1; $i++) {
            $string .= "0";
            last if ($i == $nzeros);
        }
        return($string);    
      }
    }
    if ($x!~ /1\.\d*/ && $x !~ /^0\.\d+$/ && $x !~ /^(\d)\.(\d*)e\-(\d+)$/) {die "$x\tunexpected input to sub significantdigits";}
    my $string = "0.";
    my $first;
    my $start=0;
    
    my $depth;
    if ($x !~ /e/) {
        $depth = length($x)-2
    } else {
        if ($x =~ /(\d)\.(\d*)e\-(\d+)/ || die "unexpected e format $x") {
            $depth = $3 + length($2);
        }    
    }
   
  if ($x !~ /e/) {
    for (my $i=1; $i<=$depth; $i++) {

        #print "\n$string\n";

        #my $y1 = int($x*(10**$i)); 
        ##my $y2 = int($string*(10**($i-1))); #this method sometimes gives somesort of rounding/carryover error
        ##$y2 .= "0";  
        #my $y2 = substr($string, 2)."0";
        #my $y = $y1 - $y2;

        my $y = substr($x, $i+1, 1);

        #print "$string\t";
        #print "$y1 - $y2 = $y\n";
        $string .= $y;
        #print "$string\n";

        if (!$start && $y != 0) {
            $first = $i;
            $start=1;
        }
        if ($start && $i == $first + $n - 1) {
            #my $z1 = int($x*(10**($i+1)));
            #my $z2 = substr($string, 2)."0";
            #my $z = $z1 - $z2;
            my $z = substr($x, $i+2, 1);
            #print "z $z\n";
            if ($z > 4) { #round up   #note will not add zeros beyond the last digit of the original input
                #print "$x\t$string\n";
                $string += 1/(10**$i);
                #fill with zeros if nec
                #print "$x\t$string\n";
              
                my $newfirst = 0;
                
               if ($string eq "1") {
                   $string = "1.";
               } else {  
                my $depth2;
                $depth2 = length($string)-2;

                for (my $i=1; $i<=$depth2; $i++) {
                     #my $y1 = int($string*(10**$i));

                       my $y1 = substr($string, $i+1, 1);

                     
#                     print "power ".10**$i."\n";
#                     
#                     print "mult ".$string*(10**($i+7))."\n";
#                     
#                     print "int ".int($string*(10**$i))."\n";
#                     
#
#                    my $t1 = 10**$i;
#                    my $t2 = $string * $t1;
#                    my $t3 = int(1);
#                    
#                    print "$t1 $t2 $t3\n";
#                     
                     #print "$string\t$i\t$y1\n";
                     if ($y1 != 0) {
                         $newfirst = $i;
                         last;
                     }
                }
               } 
                
                if ($newfirst != $first &&  $newfirst != $first-1) {
                    die "newfirst?\t$first\t$newfirst\t$x\t$string";
                }    
                
                my $jmax = ($newfirst + $n + 1) - length($string);
                #print "\$jmax $jmax\n";
                #print "$jmax = ($newfirst + $n + 1) - length($string)\n";
                for (my $j=1; $j<=$jmax; $j++) {
                    $string .= "0";
                    #print "$string\n";
                }
            }
            last;    
        }
    }
  } else {  
    if ($x =~ /(\d)\.(\d*)e\-(\d+)/ || die "unexpected e format $string") {
        my ($x1, $y1, $z1) = ($1, $2, $3);
        die "unexpected e format starting with a zero $string" if ($x1 == 0);
        my $nsd = length($x1) + length($y1);
        if ($nsd > $n) {
            my $a = $x1.$y1;
            my $b = substr($a, 0, $n);
            my $c = substr($a, $n, 1);
            if ($c > 4) {
                $b += 1;
            }
            if (length($b) == $n+1) {   #if rounding up added an extra digit   
                $b = substr($b, 0, $n);
                $z1 -= 1; 
            } elsif (length($b) != $n) {die;}
            $x1 = substr($b, 0, 1);
            $y1 = substr($b, 1); 
            $string = "$x1.${y1}e-$z1";
        } else {
            $string = $x;
        }
      }
  }  
    
    return ($string);

}

############################################################
#                                                          
# ttest - Student's t-test
#                     uses: betai
#                             Student t-test calculations
#
############################################################

sub ttest {

    my ($x1, $x2) = @_;

    my ($error, $mean1, $var1) = &meanvar($x1);
    die $error if $error;
    my ($error, $mean2, $var2) = &meanvar($x2);
    die $error if $error;
    
    my $n1minus1 = $#{$x1};
    my $n2minus1 = $#{$x2};
    
    #die if ($varthresh eq "");
    #if ($var1 < $varthresh) {$var1 = $varthresh;}
    #if ($var2 < $varthresh) {$var2 = $varthresh;}

    my ($score, $prob);
    if ($var1 <= 0 && $var2 <= 0) {
        die "mean1 var1 mean2 var2 $mean1 $var1 $mean2 $var2";
    } elsif ($var1 <= 0 || $var2 <= 0) {
        warn "mean1 var1 mean2 var2 $mean1 $var1 $mean2 $var2";
        $score = 0;
        $prob = 1;
    } else {
        my $n1 = $n1minus1+1;
        my $n2 = $n2minus1+1;

if (1) {
        #my $temp = ($var1/$n1minus1 + $var2/$n2minus1);
        #if ($temp<0) {
        #    print "$var1/$n1minus1 + $var2/$n2minus1";
        #}
        $score = ($mean1 - $mean2)/sqrt($var1/$n1minus1 + $var2/$n2minus1);
#golub
#        $score = ($mean1 - $mean2)/(sqrt($var1) + sqrt($var2));
} elsif (0) { #less stringent version of t-test (used by Excel and SAS)
        $score = ($mean1 - $mean2)/sqrt($var1/$n1 + $var2/$n2);
        if (0) {
            my $sd1 = sqrt($var1);
            my $sd2 = sqrt($var2);
            print "$mean1 $sd1 $var1, $mean2 $sd2 $var2);";
        }
}  else {die;} #first block above likely should be used

        my $t = $score;
        my $tails = 2;

        ##estimate of degrees of freedom
        my $df = ($var1/$n1+$var2/$n2)**2 / (($var1/$n1)**2/$n1minus1+($var2/$n2)**2/$n2minus1);

        $prob = &betai(0.5*$df, 0.5, $df/($df+$t**2));
    }
    return ($score, $prob, $mean1, $var1, $mean2, $var2);
}


############################################################
#                                                          
# ttest_w_error - Student's t-test
#                     uses: betai
#                             Student t-test calculations
#
############################################################

sub ttest_w_error {

    my ($x1, $x2) = @_;

    my $error;
    my ($mean1, $var1, $mean2, $var2);
    ($error, $mean1, $var1) = &meanvar($x1);
    return $error if $error;
    ($error, $mean2, $var2) = &meanvar($x2);
    return $error if $error;
    
    my $n1minus1 = $#{$x1};
    my $n2minus1 = $#{$x2};
    
    #die if ($varthresh eq "");
    #if ($var1 < $varthresh) {$var1 = $varthresh;}
    #if ($var2 < $varthresh) {$var2 = $varthresh;}

    my ($score, $prob);
    if ($var1 <= 0 && $var2 <= 0) {
    	$error = "mean1 var1 mean2 var2 $mean1 $var1 $mean2 $var2";
    	return($error);
    } elsif ($var1 <= 0 || $var2 <= 0) {
        warn "mean1 var1 mean2 var2 $mean1 $var1 $mean2 $var2";
        $score = 0;
        $prob = 1;
    } else {
        my $n1 = $n1minus1+1;
        my $n2 = $n2minus1+1;

if (1) {
        #my $temp = ($var1/$n1minus1 + $var2/$n2minus1);
        #if ($temp<0) {
        #    print "$var1/$n1minus1 + $var2/$n2minus1";
        #}
        $score = ($mean1 - $mean2)/sqrt($var1/$n1minus1 + $var2/$n2minus1);
#golub
#        $score = ($mean1 - $mean2)/(sqrt($var1) + sqrt($var2));
} elsif (0) { #less stringent version of t-test (used by Excel and SAS)
        $score = ($mean1 - $mean2)/sqrt($var1/$n1 + $var2/$n2);
        if (0) {
            my $sd1 = sqrt($var1);
            my $sd2 = sqrt($var2);
            print "$mean1 $sd1 $var1, $mean2 $sd2 $var2);";
        }
}  else {die;} #first block above likely should be used

        my $t = $score;
        my $tails = 2;

        ##estimate of degrees of freedom
        my $df = ($var1/$n1+$var2/$n2)**2 / (($var1/$n1)**2/$n1minus1+($var2/$n2)**2/$n2minus1);

        $prob = &betai(0.5*$df, 0.5, $df/($df+$t**2));
    }
    return ($error, $score, $prob, $mean1, $var1, $mean2, $var2);
}


############################################################
#                                                          
# betai - incomplete beta function
#                     uses: betacf, gammln
#                             Student t-test calculations
#
############################################################

sub betai {
    #returns the incomplete beta function I_x(a,b)
    #uses betacf, gammln
    my ($a, $b, $x) = @_;
    my ($betai);
    my ($bt, $betacf, $gammln);
    if ($x<0 || $x>1) {
        return("error:bad arguement x $x in sub betai");
    }
    if ($x==0 || $x==1) {
        $bt=0;
    } else {  #factors in front of the continued fraction
        $bt = exp(gammln($a+$b)-gammln($a)-gammln($b)+$a*log($x)+$b*log(1-$x));
    }
    if ($x < ($a+1)/($a+$b+2)) {   #use continued fraction directly
        $betai = $bt*betacf($a,$b,$x)/$a;
        return($betai);
    } else {
        $betai = 1-$bt*betacf($b,$a,1-$x)/$b;  #use continued fraction after making the symmetry transformation
        return($betai);
    }
}


############################################################
#                                                          
# betacf - used by betai
#
#
############################################################

sub betacf {
    #used by betai:  evaluates continued fraction for incomplete beta function by modified Lentz's method
    my ($a, $b, $x) = @_;
    my $betacf;
    my $MAXIT = 100;
    my ($EPS, $FPMIN) = (3.e-7, 1.e-30);
    my ($m, $m2);
    my ($aa, $c, $d, $del, $h, $qab, $qam, $qap);
    
    $qab = $a+$b;    #these q's will be used in factors that occur in the coefficients
    $qap = $a+1;
    $qam = $a-1;
    $c=1;
    $d=1-$qab*$x/$qap;
    if (abs($d) < $FPMIN) {$d = $FPMIN;}
    $d=1/$d;
    $h=$d;
    for ($m=1; $m<=$MAXIT; $m++) {
        $m2 = 2**$m;
        $aa = $m*($b-$m)*$x/(($qam+$m2)*($a+$m2));
        $d=1+$aa*$d;                         #One step (the even one) of the recurrence.
        if (abs($d) < $FPMIN) {$d=$FPMIN;}
        $c=1+$aa/$c;
        if (abs($c) < $FPMIN) {$c=$FPMIN;}
        $d=1/$d;
        $h=$h*$d*$c;
        $aa=-($a+$m)*($qab+$m)*$x/(($a+$m2)*($qap+$m2));
        $d=1+$aa*$d;                           #Next step of the recurrence (the odd one).
        if (abs($d) < $FPMIN) {$d=$FPMIN;}
        $c=1+$aa/$c;
        if (abs($c) < $FPMIN) {$c=$FPMIN;}
        $d=1/$d;
        $del=$d*$c;
        $h=$h*$del;
        if (abs($del-1) < $EPS) {last;}          #Are we done?
    }
    #if ($m == $MAXIT+1) {
    if ($m > $MAXIT) {
        die "a $a or b $b too big, or MAXIT $MAXIT too small in betacf";
    }
    $betacf=$h;
    return($betacf);
}


############################################################
#                                                          
# gammln - log of the gamma function
#                     uses:
#                             used by betai
#
############################################################

sub gammln {
    #Returns the value ln[gamma(xx)] for xx > 0.
    my $xx = shift;    
    die "can't evaluate log of \$xx = $xx < 0" if ($xx <=0);
    my $gammln;
    my $j;
    my ($ser,$stp,$tmp,$x,$y,@cof);
   #non-perl:Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
   #accuracy is good enough.
    @cof = (76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5);
    $stp=2.5066282746310005;
    $x=$xx;
    $y=$x;
    $tmp=$x+5.5;
    $tmp=($x+0.5)*log($tmp)-$tmp;
    $ser=1.000000000190015;
    for ($j=0; $j<=5; $j++) {
        $y++;
        $ser+=$cof[$j]/$y;
    }
    $gammln=$tmp+log($stp*$ser/$x);
    return($gammln);
}


############################################################
#                                                          
# shufflearray - randomly shuffle a vector
#
#
############################################################

sub shufflearray {
    my $x = shift;  #pointer to original vector
    my $y;          #the new shuffled vector
    my $n = $#{$x} + 1;
    my @ind;
    for (my $j=0; $j<=$n-1; $j++) {
        $ind[$j]=$j;
    }
    for (my $j=0; $j<=$n-1; $j++) {
        my $l = $n - $j;
        my $rnd = int(rand $l);
        my $newj = splice(@ind, $rnd, 1);
        $y->[$j] = $x->[$newj];
    }
    return $y;
}


############################################################
#                                                          
# arrayadd
#                                                   
############################################################

sub arrayadd {
    my ($a1, $a2) = @_;
    my $n = $#{$a1}+1;
    die if ($n != $#{$a2}+1);

    my @a3 = ();
    --$n;
    for my $i (0 .. $n) {
        $a3[$i] = $a1->[$i] + $a2->[$i];
    }
    return(\@a3);
}


############################################################
#                                                          
# arraysub
#                                                   
############################################################

sub arraysub {
    my ($a1, $a2) = @_;
    my $n = $#{$a1}+1;
    die if ($n != $#{$a2}+1);

    my @a3 = ();
    --$n;
    for my $i (0 .. $n) {
        $a3[$i] = $a1->[$i] - $a2->[$i];
    }
    return(\@a3);
}

############################################################
#                                                          
# arrayave
#                                                   
############################################################

sub arrayave {
    my ($a1, $a2) = @_;
    my $n = $#{$a1}+1;
    die if ($n != $#{$a2}+1);

    my @a3 = ();
    --$n;
    for my $i (0 .. $n) {
        $a3[$i] = ($a1->[$i] + $a2->[$i])/2;
    }
    return(\@a3);
}


############################################################
#                                                          
# arrayaveblank
#                                                   
############################################################

sub arrayaveblank {
    my ($a1, $a2) = @_;
    my $n = $#{$a1}+1;
    die if ($n != $#{$a2}+1);

    my @a3 = ();
    --$n;
    for my $i (0 .. $n) {
        if ($a1->[$i] eq "" && $a2->[$i] eq "") {
            $a3[$i] = "";
        } elsif ($a1->[$i] eq "") {
            $a3[$i] = $a2->[$i];
        } elsif ($a2->[$i] eq "") {
            $a3[$i] = $a1->[$i];
        } else {
            $a3[$i] = ($a1->[$i] + $a2->[$i])/2;
        }
    }
    return(\@a3);
}


############################################################
#                                                          
# matrixmult
#                                                   
############################################################

sub matrixmult {
    my ($m1, $m2) = @_;
    my ($rows, $cols);
    my $rows = $#{$m1}+1;
    my $cols = $#{$m2->[0]}+1;
    my @m3 = ();
    --$rows; --$cols;
    for my $i (0 .. $rows) {
        my @row = ();
        my $m1i = $m1->[$i];
        for my $j (0 .. $cols) {
            my $val = 0;
            for my $k (0 .. $cols) {
                $val += $m1i->[$k] * $m2->[$k]->[$j];
            }
            push(@row, $val);
        }
        push(@m3, \@row);
    }
    return(\@m3);
}


############################################################
#                                                          
# matrixadd
#                                                   
############################################################

sub matrixadd {
    my ($m1, $m2) = @_;
    my $rows = $#{$m1}+1;
    my $cols = $#{$m1->[0]}+1;
    die if ($rows != $#{$m2}+1);
    die if ($cols != $#{$m2->[0]}+1);

    my @m3 = ();
    --$rows; --$cols;
    for my $i (0 .. $rows) {
        for my $j (0 .. $cols) {
            $m3[$i][$j] = $m1->[$i][$j] + $m2->[$i][$j];
        }
    }
    return(\@m3);
}


#############################################################
##                                                          
## roc - receiver operator charateristic 
##       and data for graphing the roc curve
##                                                   
#############################################################
#
#sub roc { #receiver operator charateristic
#
#my ($m1, $s1, $m2, $s2) = @_;
#
#if ($m1>$m2) {
#    my $temp = $m1;
#    $m1 = $m2;
#    $m2 = $temp;
#    $temp = $s1;
#    $s1 = $s2;
#    $s2 = $temp;
#}
#
#my ($min, $max) = (1e50, -1e50);
#my @x = (($m1-4*$s1), ($m1+4*$s1), ($m2-4*$s2), ($m2+4*$s2));
#
#foreach (@x) {
#    $min = &min($_, $min);
#    $max = &max($_, $max);
#}
#
##print "min $min, max $max\n";
##warn "min ($min) < 0" if ($min < 0);
##die "max ($max) < 0" if ($max < 0);
#
#my $incr = &min($s1, $s2);
#$incr = $incr/10;
##print "incr $incr\n";
#
##OUT
##my $fileout = "roc.out";
##open(OUT, ">$fileout") || die "could not open $fileout";
#
#my $sqrt2 = sqrt(2);
#my ($lasterfc1, $lasterfc2);
#my ($area, $data);  #area = area under the x=false positives, y=true positives curve
#my $i=0;
#for (my $x=$max; $x>=$min; $x-=$incr) {
#    my $Z1 = ($x-$m1)/$s1;
#    my $erfc1 = .5 * &erfc($Z1/$sqrt2);    
#    my $Z2 = ($x-$m2)/$s2;
#    my $erfc2 = .5 * &erfc($Z2/$sqrt2); 
#    my $ratio = "-";
#    if ($erfc2 !=0) {
#        $ratio = $erfc1/$erfc2;   
#    }
#    
#    if ($lasterfc1 ne "") {
#        $area += ($erfc1 - $lasterfc1) * ($erfc2 + $lasterfc2);
#    }
#
#    ($lasterfc1, $lasterfc2) = ($erfc1, $erfc2);
#    
#    #print "$x\t$Z1\t$erfc1\t$Z2\t$erfc2\t$ratio\n";
##OUT
#    #print OUT "$erfc1\t$erfc2\n";
#    
#    $data->[$i] = [($x, $erfc1, $erfc2)];    
#    $i++;
#
#}
#$area = 0.5 * $area;  #took out of the loop for efficiency
#
##print "area  $area\n";
##OUT
##close(OUT);
#return($area, $data);
#
#}

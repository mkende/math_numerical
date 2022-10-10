package Math::Numerical;

use 5.022;
use strict;
use warnings;
use utf8;

our $VERSION = 0.02;

use feature 'signatures';
no warnings 'experimental::signatures';

use Carp;
use Config;
use English;
use Exporter 'import';
use POSIX ();
use Readonly;

=pod

=encoding utf8

=head1 NAME

Math::Numerical

=head1 SYNOPSIS

Numerical analysis and scientific computing related functions.

  use Math::Numerical ':all';

  sub f { cos($_[0]) * $_[0] ** 2 }  # cos(x)·x²

  my $root = find_root(\&f, -1, 1);

=head1 DESCRIPTION

This module offers functions to manipulate numerical functions such as root
finding (solver), derivatives, etc. Most of the functions of this module can
receive a C<$func> argument. This argument should always be a code reference (an
anonymous sub or a reference to a named code block). And that referenced
function should expect a single scalar (numeric) argument and return a single
scalar (numeric) value. For efficiency reason, it is recommended to not name the
argument of that function (see the L<example above|/SYNOPSIS>).

=head1 CONFIGURATION

For now, this module has no global configuration available. All configuration
must be directly passed to the individual functions.

=head1 FUNCTIONS

By default, none of the functions below are exported by this package. They can
be selectively imported or you can import them all with the tag C<:all>.

=cut

our @EXPORT_OK = qw(find_root bracket);
our %EXPORT_TAGS = (all => \@EXPORT_OK);

# This will need to be adapted if we start using bigfloat.
Readonly my $EPS =>
  $Config{uselongdouble} ? POSIX::LDBL_EPSILON : POSIX::DBL_EPSILON;
Readonly our $_DEFAULT_TOLERANCE => 0.00001;  # exposed for tests only.
Readonly my $DEFAULT_MAX_ITERATIONS => 100;

# Wraps the given numerical function in a way where we’re guaranteeing that it’s
# called in a scalar context and where we’re trapping its errors.
sub _wrap_func($func) {
  croak "The passed $func is not a code reference" unless ref($func) eq 'CODE';
  return sub {
    my $r = eval { &{$func} };
    return $r if defined $r;
    croak "The function failed: $EVAL_ERROR" if $EVAL_ERROR;
    croak 'The function returned no value';
  }
}

# Returns a value with the same magnitude as $val and the same sign as $sign.
sub _sign {  # ($val, $sign)
  return $_[0] * $_[1] > 0 ? $_[0] : -$_[0];
}

=head2 find_root

  find_root($func, $x1, $x2, %params)

Given a function C<$func> assumed to be continuous and a starting interval
C<[$x1, $x2]>, tries to find a root of the function (a point where the
function’s value is 0). The root found may be either inside or outside the
starting interval.

If the function is successful it returns the root found in scalar context or, in
list context, a list with the root and the value of the function at that point
(which may not be exactly C<0>).

The current implementation of this function is based on the method C<zbrent>
from the I<L<Numerical Recipes/NR>> book.

The function supports the following parameters:

=over

=item C<max_iterations>

How many iterations of our algorithm will be applied at most while trying to
find a root for the given function. This gives an order of magnitude of the
number of times that C<$func> will be evaluated. Defaults to I<100>.

=item C<do_bracket>

Whether the C<L<bracket|/bracket>> function should be used to bracket a root of
the function before finding the root. If this is set to a false value, then the
passed C<$x1> and C<$x2> values must already form a bracket (that is, the
function must take values of opposite sign at these two points). Note that, when
they do, the C<L<bracket|/bracket>> function will immediately return these
values. So, if C<find_root> return a result with C<do_bracket> set to a I<false>
value, it will return the same result when C<do_bracket> is set to a I<true>
value. However, if C<do_bracket> is set to a I<false> value, then C<find_root>
will immediately C<L<croak|Carp>> if the starting interval does not form a
bracket of a root of the function.

When set to a I<true> value, the C<L<bracket|/bracket>> function is called with
the same arguments as those given to C<find_root>, so any parameter supported by
C<L<bracket|/bracket>> can also be passed to C<find_root>.

Defaults to I<1>.

=item C<tolerance>

Defaults to I<0.00001>.

=back

In addition, as noted above, when C<do_bracket> is true, any of the parameters
supported by the C<L<bracket|/bracket>> function can be passed and they will be
forwarded to that function.

=cut

sub find_root($func, $x1, $x2, %params) {
  my $do_bracket = $params{do_bracket} // 1;
  my $tol = $params{tolerance} // $_DEFAULT_TOLERANCE;
  my $max_iter = $params{max_iterations} // $DEFAULT_MAX_ITERATIONS;
  my $f = _wrap_func($func);
  my ($xa, $xb, $xc, $xd, $xe);
  my ($fa, $fb, $fc);
  my ($p, $q, $r, $s, $tol1, $xm);
  if ($do_bracket) {
    ($xa, $xb, $fa, $fb) = bracket($func, $x1, $x2, %params);
    croak 'Can’t bracket a root of the function' unless defined $xa;
  } else {
    ($xa, $xb) = ($x1, $x2);
    ($fa, $fb) = ($f->($xa), $f->($xb));
    croak 'A root must be bracketed in [\$x1; \$x2]'
      if ($fa > 0 && $fb > 0) || ($fa < 0 && $fb <0);
  }
  ($xc, $fc) = ($xb, $fb);
  for my $i (1..$max_iter) {
    if (($fb > 0 && $fc > 0) || ($fb < 0 && $fc < 0)) {
      ($xc, $fc) = ($xa, $fa);
      $xe = $xd = $xb - $xa;
    }
    if (abs($fc) < abs($fb)) {
      ($xa, $xb, $xc) = ($xb, $xc, $xb);
      ($fa, $fb, $fc) = ($fb, $fc, $fb);
    }
    $tol1 = 2 * $EPS * abs($xb) + $tol / 2;
    $xm = ($xc - $xb) / 2;
    return wantarray ? ($xb, $fb) : $xb if abs($xm) <= $tol1 || $fb == 0;
    if (abs($xe) >= $tol1 && abs($fa) > abs($fb)) {
      $s = $fb / $fa;
      if ($xa == $xc) {
        $p = 2 * $xm *$s;
        $q = 1 - $s;
      } else {
        $q = $fa / $fc;
        $r = $fb / $fc;
        $p = $s * (2 * $xm * $q * ($q - $r)- ($xb - $xa) * ($r - 1));
        $q = ($q - 1) * ($r - 1) * ($s - 1);
      }
      $q = -$q if $p > 0;
      $p = abs($p);
      Readonly my $interp_coef => 3;
      my $min1 = $interp_coef * $xm * $q - abs($tol1 * $q);
      my $min2 = abs($xe* $q);
      if (2 * $p < ($min1 < $min2 ? $min1 : $min2)) {
        $xe = $xd;
        $xd = $p / $q;
      } else {
        $xe = $xd = $xm;
      }
    } else {
      $xe = $xd = $xm;
    }
    ($xa, $fa) = ($xb, $fb);
    if (abs($xd) > $tol1) {
      $xb +=$xd;
    } else {
      $xb += _sign($tol1, $xm);
    }
    $fb = $f->($xb);
  }
  return;
}

=head2 bracket

  bracket($func, $x1, $x2, %params)

Given a function C<$func> assumed to be continuous and a starting interval
C<[$x1, $x2]>, tries to find a pair of point C<($a, $b)> such that the function
has a root somewhere between these two points (the root is I<bracketed> by these
points). The found points will be either inside or outside the starting
interval.

If the function is successful, it returns a list of four elements with the
values C<$a> and C<$b> and then the values of function at these two points.
Otherwise it returns an empty list.

The function will C<L<croak|Carp>> if C<$x1> and C<$x2> are equal.

The current implementation of this method is a mix of the methods C<zbrac> and
C<zbrak> from the I<L<Numerical Recipes/NR>> book.

The function supports the following parameters:

=over

=item C<max_iteration>

How many iterations of our algorithm will be applied at most while trying to
bracket the given function. This gives an order of magnitude of the number of
times that C<$func> will be evaluated. Defaults to I<100>.

=item C<do_outward>

Whether the function will try to bracket a root in an interval larger than the
one given by C<[$x1, $x2]>. Defaults to I<1>.

=item C<do_inward>

Whether the function will try to bracket a root in an interval smaller than the
one given by C<[$x1, $x2]>. Defaults to I<1>.

=item C<inward_split>

Tuning parameter describing the starting number of intervals into which the
starting interval is split when looking inward for a bracket. Defaults to I<3>.

Note that the algorithm may change and this parameter may stop working or may
take a different meaning in the future.

=item C<inward_factor>

Tuning parameter describing a factor by which the inwards interval are split
at each iteration. Defaults to I<3>.

Note that the algorithm may change and this parameter may stop working or may
take a different meaning in the future.

=item C<outward_factor>

Tuning parameter describing how much the starting interval is grown at each
iteration when looking outward for a bracket. Defaults to I<1.6>.

Note that the algorithm may change and this parameter may stop working or may
take a different meaning in the future.

=back

=cut

Readonly my $DEFAULT_INWARD_SPLIT => 3;
Readonly my $DEFAULT_INWARD_FACTOR => 3;
Readonly my $DEFAULT_OUTWARD_FACTOR => 1.6;

sub _create_bracket_inward_state ($x1, $x2, $f1, %params) {
  my $s = {};
  $s->{split} = $params{inward_split} // $DEFAULT_INWARD_SPLIT;
  croak 'inward_split must be at least 2' if $s->{split} < 2;
  $s->{factor} = $params{inward_factor} // $DEFAULT_INWARD_FACTOR;
  croak 'inward_factor must be at least 2' if $s->{factor} < 2;
  @{$s}{'x1', 'x2'} = ($x1, $x2);
  $s->{f1} = $f1;
  return $s;
}

sub _do_bracket_inward ($f, $s) {
  my $dx = ($s->{x2} - $s->{x1}) / $s->{split};
  my $xa = $s->{x1};
  my ($fa, $fb) = ($s->{f1});
  for my $j (1..$s->{split}) {
    my $xb = $xa + $dx;
    $fb = $f->($xb);
    if ($fa * $fb < 0) {
      $s->{ret} = [$xa, $xb, $fa, $fb];
      return 1;
    }
    ($xa, $fa) = ($xb, $fb);
  }
  $s->{split} *= $s->{factor};
  return 0;
}

sub _create_bracket_outward_state ($f, $x1, $x2, $f1, %params) {
  my $s = {};
  $s->{factor} = $params{outward_factor} // $DEFAULT_OUTWARD_FACTOR;
  croak 'outward_factor must be larger than 1' if $s->{factor} <= 1;
  @{$s}{'x1', 'x2'} = ($x1, $x2);
  @{$s}{'f1', 'f2'} = ($f1, $f->($x2));
  return $s;
}

sub _do_bracket_outward ($f, $s) {
  if ($s->{f1} * $s->{f2} < 0) {
    $s->{ret} = [@{$s}{'x1', 'x2', 'f1', 'f2'}];
    return 1;
  }
  if (abs($s->{f1}) < abs($s->{f2})) {
    $s->{x1} += $s->{factor} * ($s->{x1} -$s->{x2});
    $s->{f1} = $f->($s->{x1});
  } else {
    $s->{x2} += $s->{factor} * ($s->{x2} -$s->{x1});
    $s->{f2} = $f->($s->{x2});
  }
  return 0;
}

sub bracket ($func, $x1, $x2, %params) {
  croak "\$x1 and \$x2 must be distinct in calls to Math::Numerical::bracket (${x1})" if $x1 == $x2;
  my $max_iter = $params{max_iterations} // $DEFAULT_MAX_ITERATIONS;
  croak 'max_iteration must be positive' if $max_iter <= 0;

  my $f = _wrap_func($func);
  my $f1 = $f->($x1);

  my $inward_state;
  if ($params{do_inward} // 1) {
    $inward_state = _create_bracket_inward_state($x1, $x2, $f1, %params);
  }
  my $outward_state;
  if ($params{do_outward} // 1) {
    $outward_state = _create_bracket_outward_state($f, $x1, $x2, $f1, %params);
  }

  croak 'One of do_outward and do_inward at least should be true'
    unless defined $outward_state || defined $inward_state;

  for my $i (1..$max_iter) {
    # We start with outward because the first iteration does nothing and just
    # checks the bounds that were given by the user.
    if (defined $outward_state && _do_bracket_outward($f, $outward_state)) {
      return @{$outward_state->{ret}};
    }
    if (defined $inward_state && _do_bracket_inward($f, $inward_state)) {
      return @{$inward_state->{ret}};
    }
    # We stop doing the inward algorithm when the number of splits exceeds
    # max_iteration, to bound the number of times the function is executed to a
    # reasonable value.
    if (defined $inward_state && $inward_state->{split} > $max_iter) {
      undef $inward_state;
    }
  }
  return;
}

=head1 AUTHOR

Mathias Kende

=head1 COPYRIGHT AND LICENSE

Copyright 2022 Mathias Kende

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=head1 SEE ALSO

=over

=item NR

L<http://numerical.recipes/>

=back

=cut

1;

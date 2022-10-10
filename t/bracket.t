use strict;
use warnings;
use utf8;

use Test2::V0;

use Readonly;
use Math::Numerical 'bracket';

use Carp;
use Test2::Tools::Compare 'validator', 'D';

$Carp::Verbose = 1;

sub float_lt {
  my ($val) = @_;
  return validator('<', $val, sub { $_ < $val });
}
sub float_gt {
  my ($val) = @_;
  return validator('>', $val, sub { $_ > $val });
}

Readonly my $PI => 4 * atan2(1, 1);

is([bracket(\&CORE::cos, 0, 1)], [float_lt($PI / 2), float_gt($PI / 2), D(), D()]);

done_testing;

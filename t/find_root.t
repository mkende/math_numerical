use Test2::V0;

use strict;
use warnings;
use utf8;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Math::Numerical 'find_root';

use Carp;
use Test2::Tools::Compare 'float';

$Carp::Verbose = 1;

use constant PI    => 4 * atan2(1, 1);

is(Math::Numerical::_DEFAULT_TOLERANCE, !number(0));

is(find_root(\&CORE::cos, 0, 3, do_bracket => 0), float(PI / 2, tolerance => 0.00001));
is(find_root(\&CORE::cos, 0, 1), float(PI / 2, tolerance => 0.00001));

done_testing;

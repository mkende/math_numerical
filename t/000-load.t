# DO NOT EDIT! This file is written by perl_setup_dist.
# If needed, you can add content at the end of the file.

#!/usr/bin/perl

use strict;
use warnings;
use Test2::V0;

BEGIN {
  ok(eval "use Math::Numerical; 1", "use Math::Numerical");
}
{
  no warnings 'once';
  note("Testing Math::Numerical $Math::Numerical::VERSION, Perl $], $^X, $ENV{SHELL}");
}

done_testing;

# End of the template. You can add custom content below this line.

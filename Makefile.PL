# DO NOT EDIT! This file is written by perl_setup_dist.
# If needed, you can add content at the end of the file.

use 5.022;
use strict;
use warnings;
use ExtUtils::MakeMaker::CPANfile;


WriteMakefile(
    NAME             => 'Math::Numerical',
    DISTNAME         => 'Math-Numerical',
    AUTHOR           => q{Mathias Kende <mathias@cpan.org>},
    VERSION_FROM     => 'lib/Math/Numerical.pm',
    ABSTRACT         => q{Algorithms for numerical functions (solver, root findings, etc.).},
    LICENSE          => 'mit',
    
    MIN_PERL_VERSION => '5.022',
    MAN3PODS         => {},
    # Directories in which we look for Makefile.PL. In general could be omitted but is needed in case there
    # is a file named Makefile.PL in a sub-directory. Should be customized if such a file needs to be
    # processed by ExtUtils.
    DIR              => [],
    NO_MYMETA => 1,
    META_MERGE => {
      'meta-spec' => { version => 2 },
      # Goes with NO_MYMETA (which would provide the dynamic config).
      dynamic_config => 0,
      no_index => {
        directory => [ 'local', 'vendor' ],
      },
      resources => {
        repository => {
          type => 'git',
          url => 'git@github.com:mkende/math_numerical.git',
          web => 'https://github.com/mkende/math_numerical',
        },
        bugtracker => {
          web => 'https://github.com/mkende/math_numerical/issues',
        },
      },
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => '.gz', },
    clean => { FILES => 'Math-Numerical-*' },
);

sub MY::postamble {
  my ($self) = @_;

  my @postamble;
  push @postamble, ::postamble() if *::postamble{CODE};

  push @postamble, <<"MAKE_FRAGMENT";
distupload: tardist
\tcpan-upload --directory Math-Numerical \$(DISTVNAME).tar\$(SUFFIX)

cover:
\tcover -test

critic:
\tperlcritic lib 

tidy:
\tperltidy -b -bext='/' \$(shell find lib -name "*.pm")

spelling: export EXTENDED_TESTING = 1
spelling:
\t\$(PERLRUN) t/001-spelling.t --interactive

alltest: export EXTENDED_TESTING = 1
alltest: test
MAKE_FRAGMENT

  return join "\n", @postamble;
}

# You can add below this template a `postamble` sub that returns more content to
# add to the generated Makefile.

# End of the template. You can add custom content below this line.

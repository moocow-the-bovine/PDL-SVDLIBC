use PDL::Core::Dev;
use ExtUtils::MakeMaker;

##-- read in user variables
do "Config.PL";
die("$0: 'do Config.PL' failed: $@") if ($@);


PDL::Core::Dev->import();
require "pdlmaker.plm";

$package = ["SVDLIBC.pd",SVDLIBC,PDL::SVDLIBC];
%hash = pdlmaker_init($package);

$hash{AUTHOR} = 'Bryan Jurish';
$hash{ABSTRACT} = 'PDL interface to Doug Rohde\'s SVD C Library';
$hash{PREREQ_PM}{PDL} = '2.4.2'; ##-- ensure PDL-ability
$hash{DIR} = [];
$hash{VERSION_FROM} = 'SVDLIBC.pd';
$hash{realclean}{FILES} .= '*~ *.tmp README.txt Config.cache las2.c svdlib.c svdutil.c svdutil.h svdlib.h';

push(@{$hash{LIBS}},    ' '.$cconfig{SVD_LIBS});
$hash{INC}           .= ' '.$cconfig{SVD_INC};
$hash{LDFLAGS}       .= ' '.$cconfig{SVD_LDFLAGS};
$hash{OPTIMIZE}      .= ' '.$cconfig{SVD_CFLAGS};
$hash{CCFLAGS}       .= ' '.$cconfig{SVD_DEFS};

$hash{OBJECT}         .= ' las2.o svdlib.o svdutil.o';


WriteMakefile(%hash);

sub MY::depend {
  package MY;
  my $inherited = shift->SUPER::depend(@_);
  $inherited .= (
		 "\n\n"
		 ##--
		 ."SVDLIBC.o: svdutil.h svdlib.h\n\n"
		 ."las2.o: las2.c svdutil.h svdlib.h\n\n"
		 ."svdlib.o: svdlib.c svdutil.h svdlib.h\n\n"
		 ."svdutil.o: svdutil.c svdutil.h svdlib.h\n\n"
		);
}

sub MY::special_targets {
  package MY;
  my $inherited = shift->SUPER::special_targets(@_);
  #$inherited =~ s/^(\.SUFFIXES\s*:.*)$/$1 .pod .rpod .pm .man .txt .html;/
  $inherited .= "

##-- SVDLIBC sources
SVDLIBC = $main::cconfig{SVD_SRCDIR}

las2.c: \$(SVDLIBC)/las2.c
\t\$(PERL) hack-int-types.perl < \$< > \$@

svdlib.c: \$(SVDLIBC)/svdlib.c
\t\$(PERL) hack-int-types.perl < \$< > \$@

svdutil.c: \$(SVDLIBC)/svdutil.c
\t\$(PERL) hack-int-types.perl < \$< > \$@

svdlib.h: \$(SVDLIBC)/svdlib.h
\t\$(PERL) hack-int-types.perl < \$< > \$@

svdutil.h: \$(SVDLIBC)/svdutil.h
\t\$(PERL) hack-int-types.perl < \$< > \$@

"
}

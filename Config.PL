## File: PDL-SVDLIBC/Config.PL
## Description: user variables for PDL::SVDLIBC package

##-- load cache
if (-e "./Config.cache") {
  do "./Config.cache";
  if ($@) {
    warn("$0: could not load cache data from 'Config.cache': $@");
  }
}

##--
## $val = cprompt($key, $message)
## $val = cprompt($key, $message, $default)
##   + sets $cconfig{$key}
sub cprompt {
  my ($key, $msg, $default)=@_;
  return $cconfig{$key} if (defined($cconfig{$key}));
  $default = '' if (!defined($default));
  print STDERR "$msg [$default] ? ";
  my $answer = <STDIN>;
  chomp($answer);
  return $cconfig{$key} = ($answer eq '' ? $default : $answer);
}

## $SVD_SRCDIR : location of SVDLIBC sources
cprompt('SVD_SRCDIR', "SVDLIBC source directory", '../SVDLIBC');

## $SVD_INC : include flags for SVDLIBC
#cprompt('SVD_INC', "include flags for SVDLIBC", '-I/usr/local/include/svdlibc');
#cprompt('SVD_INC', "Additional include flags for SVDLIBC", "-I$cconfig{SVD_SRCDIR}");
#$cconfig{SVD_INC} = "-I$cconfig{SVD_SRCDIR}";
$cconfig{SVD_INC} = "-I.";

## $SVD_DEFS : additional defines for SVDLIBC
cprompt('SVD_DEFS', "Additional defines for SVDLIBC", "-D__SVDLIBC_LONG_LONG='long long' -D__SVDLIBC_LONG='int'");

## $SVD_CFLAGS : additional compiler flags for SVDLIBC
# + do *NOT* use -ffast-math !
#cprompt('SVD_CFLAGS', "Additional compiler flags for SVDLIBC", "-O3 -pipe -march=prescott -mmmx -msse -msse2 -msse3 -mfpmath=sse");
#cprompt('SVD_CFLAGS', "Additional compiler flags for SVDLIBC", "$ENV{CFLAGS}");
#cprompt('SVD_CFLAGS', "Additional compiler flags for SVDLIBC", "-O2 -pipe");
cprompt('SVD_CFLAGS', "Additional compiler flags for SVDLIBC", $ENV{OPTIMIZE});

## $SVD_LDFLAGS : linker flags for SVDLIBC
#cprompt('SVD_LDFLAGS', "linker flags for SVDLIBC", '-L/usr/local/lib');
#cprompt('SVD_LDFLAGS', "Additional linker flags for SVDLIBC", '');

## $SVD_LIBS : libs to link
#cprompt('SVD_LIBS', "SVDLIBC libraries", '-lsvd');
#cprompt('SVD_LIBS', "Additional libraries for SVDLIBC", '');

## $SVD_OBJ : additional object files and/or static libs to link
#cprompt('SVD_OBJ', "SVDLIBC objects",
	#'/usr/local/lib/libsvd.a',
	#join(' ', qw(las2.o svdlib.o svdutil.o)
#      );


##-- save cache file
open(CACHE,">Config.cache")
  or die("$0: open failed for Config.cache: $!");
print CACHE "## Automatically generated, remove to re-configure!\n";
while (($key,$val)=each(%cconfig)) {
  print CACHE "\$cconfig{'$key'} = q($val);\n";
}
print CACHE "1;\n";
close CACHE;

1; ##-- return nicely

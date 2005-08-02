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

## $SVD_INC : include flags for SVDLIBC
cprompt('SVD_INC', "include flags for SVDLIBC", '-I/usr/local/include/svdlibc');

## $SVD_LDFLAGS : linker flags for SVDLIBC
cprompt('SVD_LDFLAGS', "linker flags for SVDLIBC", '-L/usr/local/lib');

## $SVD_LIBS : libs to link
cprompt('SVD_LIBS', "SVDLIBC libraries", '-lsvd');

## $SVD_OBJ : additional object files and/or static libs to link
cprompt('SVD_OBJ', "SVDLIBC objects",
	#'/usr/local/lib/libsvd.a',
	#join(' ', qw(las2.o svdlib.o svdutil.o)
       );

##-- save cache file
open(CACHE,">Config.cache")
  or die("$0: open failed for Config.cache: $!");
print CACHE "## Automatically generated, remove to re-configure!\n";
while (($key,$val)=each(%cconfig)) {
  print CACHE "\$cconfig{'$key'} = '$val';\n";
}
print CACHE "1;\n";
close CACHE;

1; ##-- return nicely

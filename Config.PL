## File: PDL-SVDLIBC/Config.PL
## Description: user variables for PDL::SVDLIBC package

##-- load cache
if (0 && -e "./Config.cache") {
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
 my $answer = ExtUtils::MakeMaker::prompt("<Config.PL> $msg [$default] ?");
  chomp($answer);
  return $cconfig{$key} = ($answer eq '' ? $default : $answer);
}

##--------------------------------------------------------------
## Sources

## $SVD_SRCDIR : location of SVDLIBC sources
#my $SVD_SRC0 = (grep {-d $_} qw(./SVDLIBC ../SVDLIBC))[0] // "http://tedlab.mit.edu/~dr/SVDLIBC/svdlibc.tgz";
my $SVD_SRCDIR  = 'SVDLIBC';
my $SVD_SRCFILE = 'svdlibc.tgz';
my $SVD_SRC0 = ((grep {-d $_} qw(./SVDLIBC ../SVDLIBC)),
		(grep {-f $_} ($SVD_SRCFILE, 'svdlibc.tar.gz')),
		"http://tedlab.mit.edu/~dr/SVDLIBC/svdlibc.tgz",
	       )[0];
my $SVD_SRC  = cprompt('SVD_SRC', "SVDLIBC source directory, archive, or URL", $SVD_SRC0);

$SVD_SRC =~ s{^file://}{};
if (-d $SVD_SRC) {
  ##-- local directory source
  $cconfig{SVD_SRCDIR} = $SVD_SRC;
}
elsif (-f $SVD_SRC) {
  ##-- local file source: buffer
  ;
}
elsif ($SVD_SRC =~ m{^\w+://}) {
  ##-- url source
  print STDERR "Info: attempting to download SVDLIBC source from $SVD_SRC\n";
  $SVD_SRC = get_url($SVD_SRC,$SVD_SRCFILE);
}

if (-f $SVD_SRC) {
  ##-- source archive: unpack it
  print STDERR "Info: attempting to unpack SVDLIBC source from $SVD_SRCFILE ... ";
  require Archive::Tar;
  require File::Basename;
  my $tar = Archive::Tar->new($SVD_SRCFILE)
    or die("Error: SVDLIBC source file $SVD_SRCFILE doesn't look like a tar archvie");
  $tar->contains_file('SVDLIBC/svdlib.h')
    or die("Error: SVDLIBC source archive $SVD_SRCFILE doesn't contain required file(s)");
  $tar->extract()
    or die("Error: failed to extract SVDLIBC source from $SVD_SRCFILE");
  $cconfig{SVD_SRCDIR} = './SVDLIBC';
  print STDERR "success.\n";
}

##-- sanity check
if (! -d $cconfig{SVD_SRCDIR}) {
  warn("Warning: source directory SVD_SRCDIR=$SVD_SRCDIR does not exist!");
  return;
}


##--------------------------------------------------------------
## Flags

## $SVD_INC : include flags for SVDLIBC
#cprompt('SVD_INC', "include flags for SVDLIBC", '-I/usr/local/include/svdlibc');
#cprompt('SVD_INC', "Additional include flags for SVDLIBC", "-I$cconfig{SVD_SRCDIR}");
#$cconfig{SVD_INC} = "-I$cconfig{SVD_SRCDIR}";
$cconfig{SVD_INC} = "-I.";

## $SVD_DEFS : additional defines for SVDLIBC
cprompt("ENABLE_PDL_INDX", "Enable 64-bit PDL_Indx support if available (requires PDL >= v2.007)", "y");
$cconfig{ENABLE_PDL_INDX} = ($cconfig{ENABLE_PDL_INDX} && $cconfig{ENABLE_PDL_INDX} !~ /^[n0]/ ? 1 : 0);

require PDL::Core;
my ($LONG,$LONGLONG) = $cconfig{ENABLE_PDL_INDX} && defined(&PDL::indx) ? ('PDL_Indx','PDL_Indx') : ('PDL_Long','PDL_LongLong');
cprompt('SVD_LONG',     "SVDLIBC long integer datatype", "$LONG");
cprompt('SVD_LONGLONG', "SVDLIBC long long integer datatype", "$LONGLONG");
cprompt('SVD_DEFS', "Additional defines for SVDLIBC", "-D__SVDLIBC_LONG='$cconfig{SVD_LONG}' -D__SVDLIBC_LONG_LONG='$cconfig{SVD_LONGLONG}'");

## $SVD_CFLAGS : additional compiler flags for SVDLIBC
# + do *NOT* use -ffast-math !
#my $SVD_CFLAGS = $ENV{CFLAGS};
#my $SVD_CFLAGS = "-O2 -pipe";
#my $SVD_CFLAGS = "-O3 -pipe -march=native -mtune=native -mfpmath=sse";
#cprompt('SVD_CFLAGS', "SVDLIBC compiler flags", $SVD_CFLAGS);

## $SVD_LDFLAGS : linker flags for SVDLIBC
#cprompt('SVD_LDFLAGS', "SVDLIBC linker flags", '-L/usr/local/lib');

## $SVD_LIBS : libs to link
#cprompt('SVD_LIBS', "SVDLIBC libraries", '');

## $SVD_OBJ : additional object files and/or static libs to link
#my $SVD_OBJ = '/usr/local/lib/libsvd.a';
#my $SVD_OBJ = join(' ', qw(las2.o svdlib.o svdutil.o);
#cprompt('SVD_OBJ', "SVDLIBC objects", $SVD_OBJ);

##--------------------------------------------------------------
## save cache file
use Data::Dumper;
open(CACHE,">Config.cache")
  or die("$0: open failed for Config.cache: $!");
print CACHE
  ("## Automatically generated, remove to re-configure!\n",
   Data::Dumper->Dump([\%cconfig],['*cconfig']),
   "1; ##-- be happy\n",
  );
close CACHE;

##--------------------------------------------------------------
## subs

## $outfile_or_undef = get_url($url, $outfile)
##  + URL download code adapted from HTTP::UA::Parser module's Makefile.PL
sub get_url {
  my ($url,$outfile) = @_;

  ##-- curl
  print STDERR "+ Trying curl... ";
  if (system('curl','-sS','-o',$outfile,$url)==0) {
    print STDERR "success.\n";
    return $outfile;
  }
  print STDERR "failed.\n";

  ##-- wget
  print STDERR "+ Trying wget... ";
  if (system('wget', '-nv', '-O',$outfile, $url)==0) {
    print STDERR "success.\n";
    return $outfile;
  }
  print STDERR "failed.\n";

  ##-- lwp-request
  my ($infh,$outfh);
  my $rc = 'failed';
  open($outfh,">$outfile") or die("open failed for $outfile: $!");
  print STDERR "+ Trying lwp-request... ";
  eval {
    open($infh, "-|", 'lwp-request', $url) or die("lwp-request $url failed: $!");
    print $outfh $_ while (defined($_=<$infh>));
    $rc = 'success';
  };
  print STDERR "$rc.\n";
  return $outfile if ($rc ne 'failed');

  ##-- try LWP::UserAgent
  print STDERR "+ Trying LWP::UserAgent... ";
  eval "use LWP::UserAgent";
  if ($@) {
    print STDERR
      ("failed.\n",
       "\n",
       " Couldn't locate LWP::UserAgent module\n",
       " LWP::UserAgent is required to fetch SVDLIBC sources from server\n",
       " Please install it or get SVDLIBC sources manually from\n",
       " $url\n",
       " and place it as $SVD_SRCFILE the root folder of this distro\n",
       " then run Makefile.PL again\n",
       "\n",
      );
    exit;
  }
  my $ua = LWP::UserAgent->new;
  $ua->timeout(5);
  $ua->env_proxy();
  my $response = $ua->get($url);
  if ($response->is_success) {
    print STDERR "success.\n";
    print $outfh $response->content;
    return $outfile;
  }
  else {
    die "LWP::UserAgent request aborted\n";
    exit;
  }
}


1; ##-- return nicely


#!/usr/bin/perl -w

use PDL;

print "sizeof(byte) = ", PDL::howbig(byte->[0]), "\n";
print "sizeof(short) = ", PDL::howbig(short->[0]), "\n";
print "sizeof(long) = ", PDL::howbig(long->[0]), "\n";
print "sizeof(longlong) = ", PDL::howbig(longlong->[0]), "\n";
if (defined(&PDL::indx)) {
  print "sizeof(indx) = ", PDL::howbig(indx->[0]), "\n";
}

#print "\n";
print "sizeof(float) = ", PDL::howbig(float->[0]), "\n";
print "sizeof(double) = ", PDL::howbig(double->[0]), "\n";

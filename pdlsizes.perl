#!/usr/bin/perl -w

use PDL;

print "sizeof(short) = ", PDL::howbig(short->[0]), "\n";
print "sizeof(long) = ", PDL::howbig(long->[0]), "\n";
print "sizeof(longlong) = ", PDL::howbig(longlong->[0]), "\n";
print "sizeof(float) = ", PDL::howbig(float->[0]), "\n";
print "sizeof(double) = ", PDL::howbig(double->[0]), "\n";


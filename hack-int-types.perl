#!/usr/bin/perl -w

while (<>) {
  s/\blong\b\s+long\b(?:\s+int\b)?/__SVDLIBC_LONG_LONG/g;
  s/\blong\b(?:\s+int\b)?/__SVDLIBC_LONG/g;
  print;
}

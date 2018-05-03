#!/usr/bin/env perl

use strict;
my $link=$ARGV[0];
my $file=(split/:/,$link)[1];
$file =~ /(file.+)"}/g;
my $id=$1;
print $id;

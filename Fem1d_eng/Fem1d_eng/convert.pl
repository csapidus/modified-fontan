#!/usr/bin/perl
die " I need an argument" unless $ARGV[0];
$dir=$ARGV[0];
opendir DIR, $dir;
while ( $file = readdir(DIR) )
  {
    if($file ne '.' && $file ne '..'){
      $filem=$file;
      $filem =~ tr/A-Z/a-z/;
      if ($filem ne $file) {print $file . "  " . $filem . "\n"};
      if ($filem ne $file) {rename $file, $filem};
       
    }
  }

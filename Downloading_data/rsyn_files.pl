use strict;

my $base = 'rsync --copy-links --times --verbose rsync://';
my ($url, $command);

my $out = "$ARGV[1]";

while(<>) {
   chomp;
   if (/^(ftp\..*)/) {
      $url = $1;
   }
   elsif (/:\/\/(ftp.*)/) {
      $url = $1;
   }
   else {
      print "Invalid URL: $_\n";
      next;
   }
   $command = $base . $url . " " . $out;
   print "$command\n";
   print "Downloading $url\n";
   system($command); 
}

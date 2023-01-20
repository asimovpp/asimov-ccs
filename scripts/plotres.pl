#!/usr/bin/perl -w

# Test whether gnuplot is in $PATH
$gnuplot = `which gnuplot`;
chomp($gnuplot);
if (! -e $gnuplot) {
  die("\nError: gnuplot must be in your \$PATH\n\n");
}

# Read filename of file containing residual data from command line
if($#ARGV != 0) {
  die("\nusage: plotres.pl <Residual/stdout filename>\n\n");
}
else {
  # Open file containing residual data
  $filename = $ARGV[0];
  open(INFILE,$filename)
    or die("\nError: can't open input file \"$filename\"\n\n");
}

# Read column headers
while($line = <INFILE>){
  if ($line =~ /Iter/){
    @headers = split(" ",$line);
    $nheaders = scalar(@headers);
    last;
  }
}
close(INFILE);

# Create gnuplot script
open(GNUPLOT_SCRIPT,'>', "gnuplot_script");
print GNUPLOT_SCRIPT "set logscale y\n";
print GNUPLOT_SCRIPT "set format y \"%.1e\"\n";
print GNUPLOT_SCRIPT "set xlabel \"Iteration\"\n";
print GNUPLOT_SCRIPT "set ylabel \"Residual\"\n";
print GNUPLOT_SCRIPT "plot ";
for($i=1; $i<$nheaders; $i++) {
  print GNUPLOT_SCRIPT "\"$filename\" using 1:$i+1 with lines title \"$headers[$i]\",";
}
close(GNUPLOT_SCRIPT);

# Plot with gnuplot
system("gnuplot -p gnuplot_script");

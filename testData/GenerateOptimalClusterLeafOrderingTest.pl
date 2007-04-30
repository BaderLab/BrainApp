#!/usr/bin/perl
#Generate a number of fasta formatted protein sequences with random sequence except for a
#string of A's at varying positions

$aa = "acdefghiklmnpqrstvwy";
$numSeqs = 30;
$intervalSize = 5;
$matchSize = 40;
$sortOutput=1;  #sort or random order output

%seqs;  #keep sequences in a hash to output them in random order

for($i=0;$i<$numSeqs;$i++) {
  $seqString="";
  $seqIndex=$i;
  if(length($seqIndex)<length($numSeqs)) { #pad index for sorting
    $seqIndex="0"x(length($numSeqs)-length($seqIndex)).$seqIndex;
  }
  $seqString = $seqString.">test$seqIndex\n";
  for($j=0;$j<($i*$intervalSize);$j++) {
    $rand = int(rand(length($aa)));
    $seqString = $seqString.substr($aa,$rand, 1);
  }
  for($j=0;$j<$matchSize;$j++) {
    $seqString = $seqString."A";
  }
  for($j=0;$j<(($numSeqs-$i)*$intervalSize);$j++) {
    $rand = int(rand(length($aa)));
    $seqString = $seqString.substr($aa,$rand, 1);
  }
  $seqString = $seqString."\n";
  $seqs{$seqString}=1;
}

if($sortOutput==1) {
  foreach $k (sort(keys(%seqs))) {
    print $k;
  }
}
else {
  foreach $k (keys(%seqs)) {
    print $k;
  }
}
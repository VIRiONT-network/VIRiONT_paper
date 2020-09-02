use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use File::Copy;

#Input

$reportfile = shift;
$MIXED_CUTOFF = shift;
$consensusfile = shift;
$MINCOV = shift; # Minimum coverage to not call a "N" (original : 15)

#Parameters etc

$MINCOV_OUTPUT = 10; # Minimum mean coverage to output any consensus at all (original : 25)
$MINCOV_MIXED = 20; # Minimum coverage to call mixed IUPAC bases (original : 40)
$MINCOV_INDEL = 40; # Minimum coverage to call indels (original : 40)
$INDEL_CUTOFF = 0.85; # Minimum fraction of indels to call out of frame indels.(original : 0.75)
$error = 0;

&pathogen_varcaller($reportfile, $reportfile.'_variants.txt', $consensusfile, $MINCOV, $MINCOV_MIXED, $MINCOV_INDEL, $MIXED_CUTOFF, $INDEL_CUTOFF);


sub  pathogen_varcaller {
$pileupfile = shift;
$variantfile = shift;
$consensusfile = shift;
$mincoverage = shift;
$MINCOV_MIXED = shift;
$MINCOV_INDEL = shift;
$MIXED_CUTOFF = shift;
$INDEL_CONSENSUS_CUTOFF = shift;

@symbols = ('A','T','C','G');
%IUPAC = (	A => 'A',
			T => 'T',
			G => 'G',
			C => 'C',
			U => 'U',
			AT => 'W',
			TA => 'W',
			CG => 'S',
			GC => 'S',
			AC => 'M',
			CA => 'M',
			GT => 'K',
			TG => 'K',
			AG => 'R',
			GA => 'R',
			CT => 'Y',
			TC => 'Y',
			CGT => 'B',
			AGT => 'D',
			ACT => 'H',
			ACG => 'V',
			);


open $pileupf, '<', $pileupfile or die "No input pileup file";
open $variantf, '>', $variantfile or die "No output file specified";
open $consensusf, '>', $consensusfile or die "No consensus file specified";

if (!defined($mincoverage)) {
	print "No minimum coverage specified, using default ($DEFAULT_MINCOV)\n";
	$mincoverage = $DEFAULT_MINCOV;
}


my %consensus = ();
my %votelist = ();
my @chrlist = ();
$chrold = '';
$deletion = 0;
while (defined($l = <$pileupf>)) {
	($chr, $refpos, $refbase, $cov, $votes, $quals) = split '\s+',$l;
	
	if ($chr ne $chrold) {
		$chrold = $chr;
		push @chrlist, $chr;
	}
	$refbase = uc $refbase;
	$base = $refbase;
	$votes = uc $votes;
	$votes =~ s#\^.##g;
	%indels = ();
	while ($votes =~ m#(\+)(\d+)#) {
		$indels{substr($votes,$-[0],$+[1]-$-[0] + $2 +1)}++;
		substr($votes,$-[0],$+[1]-$-[0] + $2 +1) = '';
	}
	while ($votes =~ m#(\-)(\d+)#) {
		$indels{substr($votes,$-[0],$+[1]-$-[0] + $2 +1)}++;
		substr($votes,$-[0],$+[1]-$-[0] + $2 +1) = '';
	}
	
	$totvotes = 0; 
	my %votelist = ();
	SYMBOL: foreach $symbol (@symbols) {
		if ($symbols eq $refbase) {
			next SYMBOL;
		}
		$vnum = () = ($votes =~ m#$symbol#g);
		$totvotes += $vnum;
		$votelist{$symbol} = $vnum;
	}
	
	$votelist{$refbase} = () = ($votes =~ m#[.,]#g);
	$totvotes += $votelist{$refbase};
	@baselist = sort { 
		return $votelist{$b} <=> $votelist{$a};
	} keys %votelist;
	
	if ($totvotes > 0) {	
		$ratio = $votelist{$baselist[1]} / $totvotes;
	} else {
		$cov = 0;
		$ratio = 0;
		$totvotes = 1;
	}	
	if ($cov < $mincoverage) {
		$base = 'N';	
	} elsif ($ratio < $MIXED_CUTOFF || $cov < $MINCOV_MIXED) {
		$base = $baselist[0];		
	} else {
		$base = $IUPAC{$baselist[0].$baselist[1]};	
	}
	
	# SNP
	#if ($base ne $refbase && $base ne 'N') {
	#if ($base =~ m#[ATGC]#) {
		#$ratio2 =  $votelist{$baselist[0]} / $totvotes;	
	#} else {
		$ratio2 =  $votelist{$baselist[1]} / $totvotes;
	#}
	print {$variantf} $chr.' '.$refpos.' '.'cov: '.$cov.' '.$refbase.'/'.$base;
	foreach $b (@baselist) {
		print {$variantf} " $b:".$votelist{$b};
	}
	$pr = sprintf "%.1f",$ratio2 * 100;
	print {$variantf} " = $pr%\n";
	#}
	# INDEL
	$insert = '';
	KEYS: foreach $k (keys %indels) {
		if ($k =~ m#(\d+)#) {
			$indlen = $1;
		} else {
			next KEYS;
		} 
		if ($indlen < 1) {
			$indlen = 1;
		}
		
		@cutofflist = (1,1,0.8,0.3,0.7,0.6,0.3);
		if ($indlen < 7) {
			$cutoff = $INDEL_CONSENSUS_CUTOFF * $cutofflist[$indlen];
		} elsif ($indlen % 3 == 0) {
			$cutoff = $INDEL_CONSENSUS_CUTOFF * 0.3;
		} else {
			$cutoff = $INDEL_CONSENSUS_CUTOFF * 0.6;
		}
		# $cutoff = $INDEL_CONSENSUS_CUTOFF * (1/$indlen);
		
		if ($cutoff < $MIXED_CUTOFF) {
			$cutoff = $MIXED_CUTOFF;
		}
		
		if ($indels{$k} > $totvotes * $cutoff && $cov >= $MINCOV_INDEL) {
			$pr = sprintf "%.1f",100 * ($indels{$k}/$totvotes);
			print {$variantf} $chr.' '.$refpos.' '.'cov: '.$cov.' '.$refbase.' INDEL: '.$k." $pr%\n";
			if ($k =~ m#^\+\d+(.*)#) {
				$insert = $1;
			} elsif ($k =~ m#^\-(\d+).*#) {
				$deletion = $1;
			}
			
			
		}
		
	}
	
	# INDELS IN CONSENSUS
	if ($deletion > 0) {
		print "Deletion $base in consensus of $chr at $refpos\n";
		$deletion--;	
	} else {
		$consensus{$chr} = $consensus{$chr}.$base;
	} 
	if ($insert) {
		print "Insert $insert in consensus of $chr at $refpos\n";
		$consensus{$chr} = $consensus{$chr}.$insert;
	}
}

foreach $ch (@chrlist) {
	$out1 = '>'.$ch."\n";
	$out2 = uc($consensus{$ch})."\n";
	print {$consensusf} $out1;
	print {$consensusf} $out2;
}

print {$variantf} 'EOF@EOF';
close $variantf;
print "Wrote variant file: $variantfile\n";
close $consensusf;
print "Wrote consensus file: $consensusfile\n";
close $pileupf;
	
}



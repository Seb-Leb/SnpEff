#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Parse and modify NCBI's GFF (RefSeq) to make is suitable for building
# a SnpEff database
#
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

use strict;

# Debug mode?
my($debug) = 0;

my(@keys) = ();		# GFF attribute key order

#-------------------------------------------------------------------------------
# Parse GFF line
#-------------------------------------------------------------------------------
sub parseGffLine($) {
	my($l) = @_;
	my($chr, $source, $type, $start, $end, $score, $strand, $phase, $attrs);

	chomp $l;
	($chr, $source, $type, $start, $end, $score, $strand, $phase, $attrs) = split /\t/, $l;
	print STDERR "Parsing GFF line:\t$chr\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$attrs\n" if $debug;

	#---
	# Translate types
	#---
	if( $type eq 'primary_transcript' ) { $type = 'transcript'; }

	#---
	# Parse attributes as key-value pairs
	#---
	my($key, $value, $kv, %attr, @kvs);
	my(%attr) = ();
	my(@kvs) = split /\;/, $attrs;
	@keys = ();
	foreach $kv ( @kvs ) {
		($key, $value) = split /=/, $kv;
		$attr{$key} = $value;
		push( @keys, $key );
		print STDERR "\tattr{$key} = $value\n" if $debug;
	}

	return ($chr, $source, $type, $start, $end, $score, $strand, $phase, \%attr);
}

#-------------------------------------------------------------------------------
# Read chromosome ID to chromosome name mapping from file
#-------------------------------------------------------------------------------
sub readChrMap($) {
	my($chrMapFile) = @_;
	my($l, $chrName, $chrId, %chrMap);
	open CHRMAP, $chrMapFile;
	while( $l = <CHRMAP> ) {
		chomp $l;
		($chrName, $chrId) = split /\t/, $l;
		$chrMap{$chrId} = $chrName;
	}
	close CHRMAP;
	return %chrMap;
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Parse command line arguments
die "Usage: $0 file.gff id_map.txt id_map_protein.txt chr_map.txt\n" if $#ARGV < 0;
my($gff) = $ARGV[0];
die "Missing command line argument 'file.gff'\n" if $gff eq '';
my($idMapFile) = $ARGV[1];
die "Missing command line argument 'id_map.txt' (output file)\n" if $idMapFile eq '';
my($idProtMapFile) = $ARGV[2];
die "Missing command line argument 'id_map_protein.txt' (output file)\n" if $idProtMapFile eq '';
my($chrMapFile) = $ARGV[3];
die "Missing command line argument 'chr_map.txt' (input file)\n" if $chrMapFile eq '';

my(%chrMap) = readChrMap($chrMapFile);

#---
# Parse GFF files from STDIN
#---

my($idsToChange) = {
	'gene' => 1
	,'C_gene_segment' => 1
	,'D_gene_segment' => 1
	,'J_gene_segment' => 1
	,'V_gene_segment' => 1
	, 'lnc_RNA' => 1
	, 'mRNA' => 1
	, 'ncRNA' => 1
	, 'primary_transcript' => 1
	, 'rRNA' => 1
	, 'tRNA' => 1
	, 'transcript' => 1
};

my($chr, $chrid, $source, $type, $start, $end, $score, $strand, $phase, $attr);	# GFF parsed values
my(%name2id);	# Map name to IDs
my(%idOld2New);	# Map old ID to new ID
my(%prot);		# Map protein ID to transcript ID
my($l, $id, $key, $value, $name, $paren, $newParent);

print STDERR "Parsing GFF file '$gff'\n";
if( $gff =~ /.gz$/ ) {
	open GFF, "gunzip -c $gff |" || die "Cannot open file '$gff'\n";
} else {
	open GFF, $gff || die "Cannot open file '$gff'\n";
}

while( $l = <GFF> ) {
	# Skip headers
	if( $l !~ /^#/ ) {
		($chr, $source, $type, $start, $end, $score, $strand, $phase, $attr) = parseGffLine($l);

		# Map chromosme id to name
		if( $chrMap{$chr} ne '' ) {
			$chr = $chrMap{$chr};
		}

		$name = $attr->{'Name'};
		$id = $attr->{'ID'};
		$paren = $attr->{'Parent'};

		if( $idsToChange->{$type} ) {
			if( $name ne '' ) {
				# Name to ID map and create new unique ID using 'name'
				if( $name2id{$name} ne '' ) {
					$name2id{$name} .= "\t$id";
					$idOld2New{$id} = "$name.$id"; # Append 'id' to make it unique
            
					print STDERR "Duplicated name '$name' using '$idOld2New{$id}'\n" if $debug;
				} else {
					$idOld2New{$id} = $name;
					$name2id{$name} = $id;
				}
			} else {
				# No 'name' => Use same id
				$idOld2New{$id} = $id;
			}


			# Replace ID
			$id = $attr->{'ID'} = $idOld2New{$id};
		}

		#---
		# Replace 'Parent'
		#---
		if( $paren ne '' ) {
			my($oldParent, $newParent);
			$oldParent = $paren;
			$newParent = $idOld2New{$oldParent};
			if( $newParent ne '' )	{ $attr->{'Parent'} = $newParent; }
			else					{ print STDERR "WARNING: Cannot find ID for parent '$oldParent' in line: $l\n"; }

			#---
			# Protein - Transcript map
			#---
			if( $type eq 'CDS' ) {
				my($protId) = $name;
				my($trId) = $paren;
				if( $protId ne '' ) {
					$prot{$trId} = $protId;
					$prot{$newParent} = $protId;
				}
			}
		}

		# Rebuild attributes
		my($attrs) = "";
		foreach $key ( @keys ) { 
			$attrs .= ";" if $attrs ne ''; 
			$attrs .= "$key=" . $attr->{$key};
		}

		# Show GFF line
		print "$chr\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\t$attrs\n";
	}
}

close GFF;

#---
# Write IDs map file
#---
print STDERR "Creating IDs transcript map file '$idMapFile'\n";
open IDS, "> $idMapFile";
foreach $key ( sort keys %idOld2New ) { print IDS "$key\t$idOld2New{$key}\n"; }
close IDS;

print STDERR "Creating protein IDs map file '$idProtMapFile'\n";
open IDS, "> $idProtMapFile";
foreach $key ( sort keys %prot ) { print IDS "$key\t$prot{$key}\n"; }
close IDS;


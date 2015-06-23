#!/usr/local/bin/perl -w
# ==============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#             National Institute of Environmental Health Sciences
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Institute of Environmental Health 
#  Sciences (NIEHS) and the U.S. Government have not placed any restriction
#  on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NIEHS and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NIEHS and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
#===============================================================================
#
#         FILE:  mouseGT_SeqPrep.pl
#
#        USAGE:  ./mouseGT_SeqPrep.pl  -i fin -s fsnp -o fot -d dir
#
#  DESCRIPTION:  Convert user input into GenBank format;
#                Extract SNP sequence from dbSNP by rsID;
#                BLAT SNP seq to user seq, and record SNP as seq feature.
#
#      OPTIONS:
#                -i input seq file in Simple Seq Format
#                -s snp file name
#                -o output file name
#                -d temp directory
# REQUIREMENTS:  Bioperl
#         BUGS:  ---
#        NOTES:  merge seqGB_prep.pl & seqGB_snpLocate.pl
#       AUTHOR:  Mr. Hong Xu (hxu), hxu DOT hong AT gmail DOT com
#      COMPANY:  NIEHS
#      VERSION:  1.2
#      CREATED:  04/10/2010 9:23:30 AM
#     REVISION:  
#                1.1  (04/19/2012)
#                     @ add primer design input in simple seq format
#                     @ move mutation type from seqobj->desc to seqobj->accession_number
#                     @ put primer design input into seqobi->desc
#                1.2  (04/23/2012)
#                     @ generate design parameters based on sequence input
#                     
#===============================================================================

use strict;
use lib './modules';
use Getopt::Std;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use HTTP::Request::Common qw(POST GET);
use LWP::UserAgent;
use HTML::FormatText;

# get command line arguments
my %args;
getopt("isod", \%args);
if (!$args{i}) {
  die "Please specify input file name!";
}
else {
  open(FIN, "<$args{i}") or die "Cannot open input file - $args{i}";
}
if (!$args{s}) {
  die "Please specify SNP list file name!";
}
else {
  open(FSP, "<$args{s}") or die "Cannot open SNP list file - $args{s}";
}
if (!$args{o}) {
  die "Please specify output file name!";
}
my $seqio = Bio::SeqIO->new(-file => ">$args{o}",
                            -format => 'GenBank');

my $dtop = $args{d};
my $ftmp = $dtop . '/' . md5_hex(rand().$$.time()); # set temp fasta file name
my $gbio = Bio::SeqIO->new(-file => ">$ftmp.gb",
                           -format => 'GenBank');
open(FLS, ">$ftmp.log1") or die "Cannot create log file - $ftmp.log1";


# read input seq in Simple Seq Format, write in GenBank format
my $ct = 0; #DEBUG
while (my $ln = <FIN>) {
  # chomp($ln); # (hxu 06/20/12: comment out for line remove for both UNIX or WINDOWS type)
  $ln = lnClear($ln); # (hxu 06/20/12: comment out for line remove for both UNIX or WINDOWS type)
  next unless $ln; # ignore blank line
  my @aln = split(/\t/, $ln);

  # if first is entry name, initialize seq object
  if ($aln[0] eq 'Entry.name') {
    my $seq_obj = Bio::Seq->new(-display_id => $aln[1] . '|' . $aln[2]);
    my $seq_str = '';
    my $seq_cnt = 1;
    while (my $iln = <FIN>) {
      # chomp($iln); # (hxu 06/20/12: comment out for line remove for both UNIX or WINDOWS type)
      $iln = lnClear($iln); # (hxu 06/20/12: comment out for line remove for both UNIX or WINDOWS type)
      next unless $iln; # ignore blank line
      my @ailn = split(/\t/, $iln);
      
      # if entry type
      $seq_obj->accession_number($ailn[1]) if ($ailn[0] eq 'Entry.type');
      
      # if sequence feature
      if ($ailn[0] =~ m/^Sequence/) {
        $ailn[0] =~ s/^Sequence\.//;
        my $slen = length($ailn[1]);
        my $seq_nd = $seq_cnt + $slen - 1;
        my $sfeat = Bio::SeqFeature::Generic->new(-start => $seq_cnt,
                                                  -end => $seq_nd,
                                                  -primary => 'Fragment',
                                                  -tag => { name => $ailn[0] });
        $seq_str .= $ailn[1];
        $seq_obj->add_SeqFeature($sfeat);
        $seq_cnt += $slen;
      }
      
      # if entry end
      if ($ailn[0] eq 'Entry.end') {
        $seq_obj->seq($seq_str);
        $gbio->write_seq($seq_obj);
        $ct++; #DEBUG
        print FLS ("READ sequence number $ct ...\n"); #DEBUG
        last;
      }
    }
  }
}
$gbio->flush(); # flush  file handle
$gbio->DESTROY(); # close file handle

# convert GenBank sequence to Fasta format
my %feats;
my $gbt1 = Bio::SeqIO->new(-file => "$ftmp.gb",
                           -format => 'GenBank');
my $faio = Bio::SeqIO->new(-file => ">$ftmp.fa",
                           -format => 'fasta');

while (my $seq = $gbt1->next_seq) {
  $faio->write_seq($seq);
  #initialize feature 
  my @atmp = ();
  my $kid = $seq->primary_id();
  $feats{$kid} = \@atmp;
}
$gbt1->flush(); # flush file handle
$gbt1->DESTROY(); # close file handle
$faio->flush(); # flush file handle
$faio->DESTROY(); # close file handle

# get SNP info by rsSNP ID
open(FTS, ">$ftmp.snp") or die "Cannot create SNP file - $ftmp.snp";

# parse SNP file & put the result in an array
my @records = ();
while (my $ln = <FSP>) {
  $ln = trim($ln); # remove trail space
  next unless $ln;
  push @records, $ln;
}

# get SNP from NCBI
# and convert the SNP information into format for BLAT search
my $ua = LWP::UserAgent->new(); # get web agent
my %snpcnt = ();
while (scalar(@records) > 0) {
  my $rcd =shift(@records);

  my $pt = ' ';
  my $rsid;
  if ($rcd =~ /^rs(\d+)$/) {
    my $rsid = $1;
    # if rsSNP get from dbSNP
    my $req = POST 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
                [ 'db' => 'snp', 'report' => 'FASTA', 'id' => $rsid ];
    my $html = $ua->request($req)->content;
    # remove the HTML markup
    # $pt = HTML::FormatText->format_string($html);
    $pt = $html;
    # print $pt, "\n"; #DEBUG
  }
  else {
    print FLS join("\t", $rcd, $pt), "\n";
    sleep 1;
    next;
  }

  $pt = trim($pt);
  # if the sequence cannot get by network problem, put back into array and try again
  if ($pt =~ /^500/) {
    push @records, $rcd;
    sleep 1;
    next;
  }
  # temporary problem: cannot get the document
  elsif ($pt =~ /ERROR: Can't retrieve docsum/) {
    if ( (exists $snpcnt{$rcd}) && ($snpcnt{$rcd} > 5) ) {
      $pt =~ s/\n+/ /g unless ($pt eq ' ');
      print FLS join("\t", $rcd, $pt), "\n";
      sleep 1;
      next;
    }
    else {
      $snpcnt{$rcd}++;      
      push @records, $rcd;
      sleep 1;
      next;
    }  
  }
  # sequence not available from database, put into another list
  elsif ($pt =~ /^SEQUENCE CURRENTLY UNAVAILABLE/) {
    # replace line with space
    $pt =~ s/\n+/ /g unless ($pt eq ' ');
    print FLS join("\t", $rcd, $pt), "\n";
    sleep 1;
    next;
  }
  # (hxu 06/21/2012): NCBI send "empty" record like:
  #     >gnl|dbSNP|rs46290358 rs=46290358|pos=1|len=1|taxid=10090|mol=""|class=|alleles=""|build=
  #
  #     N
  elsif ($pt =~ /build=\D+/) {
    # replace line with space
    $pt =~ s/\n+/ /g unless ($pt eq ' ');
    print FLS join("\t", $rcd, $pt), "\n";
    sleep 1;
    next;
  }
  # else ...
  else {
    # $pt =~ s/1:.+//; # remove the 1st line of the FASTA file

    $pt = trim($pt); # remove leading & trailing white spaces
    $pt = trimN($pt); # remove the leading & trailing N

    # reformat DNA for BLAT search
    $pt =~ /pos=(\d+)/;
    my $ipos = $1;
    $pt =~ /len=(\d+)/;
    my $ilen = $1;
    $pt =~ /alleles="(.+)"/;
    my $allele = $1;
    $pt =~ s/^>.+//; # remove the info line of the FASTA file
    $pt =~ s/\s+//g; # remove remain white spaces

    # get 100 bp flanking sequence on both sides of the SNP
    my $Lft_st = 0;
    my $Lft_lg = $ipos - 1;
    if ($ipos > 100) {
      $Lft_st = $Lft_lg - 100;
      $Lft_lg = 100;
    } 
    my $Lft_seq = substr($pt, $Lft_st, $Lft_lg);
    my $Rgt_lg = $ilen - $ipos;
    $Rgt_lg = 100 if ($Rgt_lg > 100);
    my $Rgt_seq = substr($pt, $ipos, $Rgt_lg);

    # print left & right flanking sequences of the SNP
    print FTS ">$rcd,LFT,$allele\n$Lft_seq\n",
              ">$rcd,RGT,$allele\n$Rgt_seq\n\n";
           
    # sleep 2 sec for network courtesy 
    sleep 2;
  }
}

# do the BLAT search
my $noHead = '-noHead';
my $outfmt = '-out=pslx';
my $btmp = $ftmp . '.blat';
system("./bin/blat", "$ftmp.fa", "$ftmp.snp", $noHead, $outfmt, $btmp);
print FLS "Run BLAT ......\n"; #DEBUG blat

# parse BLAT data
# put needed data into array
my @unsort;
open(FBLT, "<$btmp") or die "Cannot open BLAT temp file - $btmp!";
while(<FBLT>) {
  my @arow = split(/\s+/, $_);
  my @sid  = split(',', $arow[9]);
  next if ($arow[17] > 1); # skip if the alignment split

  if ($arow[15] > $arow[16]) {
    ($arow[15], $arow[16]) = ($arow[16], $arow[15]);
  }
  my $score = ( $arow[0] - $arow[1] - $arow[4] * $arow[5] - $arow[6] * $arow[7] ) / ($arow[10] - $arow[3]);
  if ($score >= 0.50) {
    push @unsort, [ @sid, @arow[13,15,16,8], $score ];
  }
  # print FLS "Parse BLAT result ......\n"; #DEBUG blat
}

# data sorting 
# PRESUME that alignment be good (only one location for each SNP)
my @sorted = sort { $a->[3] cmp $b->[3] || $a->[0] cmp $b->[0] || $a->[4] <=> $b->[4] } @unsort;
my %hsnp;
for (my $i = 0; $i < scalar(@sorted); $i++) {
  unless (exists $hsnp{$sorted[$i]->[0]}) {
    $hsnp{$sorted[$i]->[0]} = $sorted[$i];
  }
  else {
    # calculate SNP location
    my $loc_prv = $hsnp{$sorted[$i]->[0]};
    my $loc_crt = $sorted[$i];
    if ($loc_prv->[1] ne $loc_crt->[1]) {
      if ( $loc_prv->[5] + 1 == $loc_crt->[4] ) {
        my $sfeat = Bio::SeqFeature::Generic->new(-start => $loc_crt->[4],
                                                  -end => $loc_crt->[4],
                                                  -primary => 'SNP',
                                                  -tag => { name => $loc_crt->[0], allele => $loc_crt->[2]});
        push @{$feats{$loc_crt->[3]}}, $sfeat;
      }        
    }
    else {
      $hsnp{$sorted[$i]->[0]} = $sorted[$i];
    }    
  }
}

# write new seq features into new GenBank file
my $gbt2 = Bio::SeqIO->new(-file => "$ftmp.gb",
                           -format => 'GenBank');
while (my $seq2 = $gbt2->next_seq) {
  # add SNP features
  my $kid = $seq2->display_id();
  my @afs = @{$feats{$kid}} if (defined @{$feats{$kid}});
  foreach my $sf ( sort {$a->start() <=> $b->start()} @afs ) {
    $seq2->add_SeqFeature($sf);
  }
  primer3Parameter($seq2); ##### generate primer3 parameters
  $seqio->write_seq($seq2);
  print FLS "Write GenBank seq in file ...\n"; #DEBUG
}
$gbt2->flush(); # flush file handle
$gbt2->DESTROY(); # close file handle
$seqio->flush(); # flush file handle
$seqio->DESTROY(); # close file handle

# remove temporary files
# system("rm $ftmp.*"); #DEBUG


###############################################
# SUB ROUTINE
###############################################
# (hxu 06/20/12: clear end of line)
sub lnClear {
  my @out = @_;
    for (@out) {
        s/[\r\n]+$//mg;          # remove end line break
    }
    return wantarray
        ? @out # many to return
        : $out[0];   # or only one
}
 

sub trim {
    my @out = @_;
    for (@out) {
        s/^\s+//mg;          # trim left
        s/\s+$//mg;          # trim right
    }
    return wantarray 
        ? @out # many to return
        : $out[0];   # or only one
}

sub trimN {
    my @out = @_;
    for (@out) {
        s/^N+//i;          # trim N left
        s/N+$//i;          # trim N right
    }
    return wantarray 
        ? @out # many to return
        : $out[0];   # or only one
}

sub primer3Parameter {
  my $pseqobj = shift;

  my $ppid = $pseqobj->display_id();
  my ($psgen, $pstyp) = split('\|', $ppid);
  $pstyp = lc($pstyp);
  my $pmutype = lc($pseqobj->accession_number);

  my @psfs = $pseqobj->get_SeqFeatures();
  # find left, right common sequence, or point mutation features
  my ($LCfeat, $PTmut, $RCfeat);
  my @XCfeats;
  foreach my $psf (@psfs) {
    my @tgvals = $psf->get_tag_values('name');
    foreach my $tgv (@tgvals) {
      $tgv = lc($tgv);
      if ($tgv eq 'left') {
        $LCfeat = $psf;
      }
      elsif ($tgv eq 'right') {
        $RCfeat = $psf;
      }
      elsif ($tgv eq 'allele') {
        $PTmut = $psf;
      }
      elsif ($tgv eq 'exclude') {
        push @XCfeats, $psf;
      }
    }
  }

  ########################### Design primer ###########################
  # mouse mutations:
  #   1. partial gene "deletion".
  #   2. partial gene "flox" & cre-induced mutation
  #   3. xeno-gene "insertion"
  #   4. gene "point" mutation
  #   5. partial gene "replace"
  #   
  # handle different mutation types differently
  #
  # Thoughts:
  #   1. deletion: design primer overlap JUNCTION of mutant first,
  #      then design primer around TARGET in wild type by using common primer from mutant primers
  #   2. flox: design primer using the 3-primer strategy of Leneuve et al. 
  #      (BioTechniques 2001)
  #   3. insertion: design primer overlap JUNCTION of mutant first,
  #      then design primer around TARGET in wild type by using common primer from mutant primers
  #   4. point: design primer ends (left or right) at POINT in wild type or mutant first,
  #      then choose the best 3 pairs to add Tm-shift tail at 5'-end.
  #   5. replace: design primer around TARGET of mutant first,
  #      then design primer around TARGET in wild type by using common primer from primer pair 1
  #####################################################################
  
  my $pfdesc;
  if ($pstyp eq 'mutant') {
    $pfdesc = 'ORDER=0; ';
    # (hxu 08/15/12: move the setup of primer targeting region here)
    # (hxu 08/07/12: added for flox-cre condition knockout mice
    #                the mutant should have the both left and right flanking sequence)
    if ( $pmutype eq 'flox' ) {
      $pfdesc = 'ORDER=1; '; # (hxu 09/05/12: for flox mice, design wild type primer first)
      unless ( $LCfeat && $RCfeat ) {
        die "Need both left and right flanking sequences";
      }
      my $target = 'TARGET=';
      $target .= ($LCfeat->location->end - 1) . ',34';
      $target .= ' ' . ($RCfeat->location->start - 1) . ',34;';
      $pfdesc .= $target;
    }
  }
  else {
    $pfdesc = 'ORDER=1; ';
    # (hxu 08/15/12: move the setup of primer targeting region here)
    # (hxu 08/07/12: added for flox-cre condition knockout mice)
    #                the mutant should have the both left and right flanking sequence
    if ( $pmutype eq 'flox' ) {
      $pfdesc = 'ORDER=0; '; # (hxu 09/05/12: for flox mice, design wild type primer first)
      unless ( $LCfeat && $RCfeat ) {
        die "Need both left and right flanking sequences";
      }
      my $target = 'TARGET=';
      $target .= ($LCfeat->location->end - 1) . ',2';
      $target .= ' ' . ($RCfeat->location->start - 1) . ',2;';
      $pfdesc .= $target;
    }
  }  
  if ($pmutype eq 'point') {
    my $pst = $PTmut->location->start;
    my $pnd = $PTmut->location->end;
    $pfdesc .= 'POINT=' . $pst . '-' . $pnd;
  }
  elsif ( ($pmutype eq 'insertion') || 
          ($pmutype eq 'deletion') ||
          ($pmutype eq 'replace') ) {
    my $target = 'TARGET=';
    if ($LCfeat) {
      $target .= ($LCfeat->location->end - 2) . ',4';
      if ($RCfeat) {
        $target .= ' ' . ($RCfeat->location->start - 2) . ',4';
      }
    }
    elsif ($RCfeat) {
      $target .= ($RCfeat->location->start - 2) . ',4';
    }
    else {
      die "No common sequence";
    }
    $target .= ';';
    $pfdesc .= $target;
  }
  elsif ( $pmutype eq 'flox' ) {
    # (hxu 08/15/12: placeholder of "if logic control" for skipping 'flox' mutation type)
  }
  else {
    die "mutation type not known";
  }
  
  # exclude
  my $excld = ' EXCLUDED=';
  for (my $idx = 0; $idx < scalar(@XCfeats); $idx++) {
    my $ft = $XCfeats[$idx];
    my $ftst = $ft->location->start;
    my $ftln = $ft->location->end - $ft->location->start;
    unless ($idx) {
      $excld .= $ftst . ',' . $ftln;
    }
    else {
      $excld .= ' ' . $ftst . ',' . $ftln;
    }
  }
  
  if (scalar(@XCfeats) > 0) {
    $pfdesc .= $excld;
  }
  else {
    # remove last ';'
    chop($pfdesc);
  }
  
  $pseqobj->desc($pfdesc);
}


=head1 Input Files

=head2 Input seq in simple sequence format (SSF)

Entry.name	Lyz2-Cre	wild
Entry.type	replace

Sequence.left	tatttcacagcagcattgcagactagctaaaggcagaagggagagactctggagtgtctcaaatgtttaagcattttaaaataaatactgtatgcttattCTTGGGCTGCCAGAATTTCTCTCATCACATAAATGAAGAAGGAAGATCAAGTGCTGAAGTCCATAGATCGGTAGGAACTTCCTGTTTTGCACACAGCTCAAATGTAGGAAACCACAAGCTGTTGGGAAAGGAGGGACTTGGAGGATGCTTAAATAGCAGGCATGCTTTCTCTAGTCAGCCAGCAGCTGACCCAGCCTCCAGTCACCATG
Sequence.change	aagactctcctgactctgggactcctcctgctttctgtcactgctcaggccaaggtctatgaacgttgtgagtttgccagaactctgaaaaggaatggaatggctggctactatggagtcagcctggccgactgtaagtctcttcttgatggctccagctaatcgctcagaacagggacaggggagcagtgagtgaccagagaggaaactgcagcatacactgtcatggttgttccctgtgtgactgaggtcatttgacagaaagcatcagtttttgtcttttgaattagagatgctggggtgagggcttggaagggacattgtgtgccatggtaaggcaacaggacctgaggggctcagggactgtcactaaattttctgttgcttcagactcaacatgagtaagggggctgagaagagatgggagagtcagccagagatcctggactgacatggtatttgtcttcacttctcttaccgtggtgtgtaacttttgagat

Entry.end


Entry.name	Lyz2-Cre	mutant
Entry.type	replace

Sequence.left	tatttcacagcagcattgcagactagctaaaggcagaagggagagactctggagtgtctcaaatgtttaagcattttaaaataaatactgtatgcttattCTTGGGCTGCCAGAATTTCTCTCATCACATAAATGAAGAAGGAAGATCAAGTGCTGAAGTCCATAGATCGGTAGGAACTTCCTGTTTTGCACACAGCTCAAATGTAGGAAACCACAAGCTGTTGGGAAAGGAGGGACTTGGAGGATGCTTAAATAGCAGGCATGCTTTCTCTAGTCAGCCAGCAGCTGACCCAGCCTCCAGTCACCATG
Sequence.change	cccaagaagaagaggaaggtgtccaatttactgaccgtacaccaaaatttgcctgcattaccggtcgatgcaacgagtgatgaggttcgcaagaacctgatggacatgttcagggatcgccaggcgttttctgagcatacctggaaaatgcttctgtccgtttgccggtcgtgggcggcatggtgcaagttgaataaccggaaatggtttcccgcagaacctgaagatgttcgcgattatcttctatatcttcaggcgcgcggtctggcagtaaaaactatccagcaacatttgggccagctaaacatgcttcatcgtcggtccgggctgccacgaccaagtgacagcaatgctgtttcactggttatgcggcggatccgaaaagaaaacgttgatgccggtgaacgtgcaaaacag

Entry.end

=head2 Input SNPs

rs37469317
rs37015784
rs36753145
rs36436921
rs38759248

=cut


=head1 Intermediate Files

=head2 Intermediate seqs in GenBank format

LOCUS       Lyz2-Cre|wild          809 bp    dna     linear   UNK 
ACCESSION   replace
FEATURES             Location/Qualifiers
     Fragment        1..309
                     /name="left"
     Fragment        310..809
                     /name="change"
ORIGIN      
        1 tatttcacag cagcattgca gactagctaa aggcagaagg gagagactct ggagtgtctc
       61 aaatgtttaa gcattttaaa ataaatactg tatgcttatt cttgggctgc cagaatttct
      121 ctcatcacat aaatgaagaa ggaagatcaa gtgctgaagt ccatagatcg gtaggaactt
      181 cctgttttgc acacagctca aatgtaggaa accacaagct gttgggaaag gagggacttg
      241 gaggatgctt aaatagcagg catgctttct ctagtcagcc agcagctgac ccagcctcca
      301 gtcaccatga agactctcct gactctggga ctcctcctgc tttctgtcac tgctcaggcc
      361 aaggtctatg aacgttgtga gtttgccaga actctgaaaa ggaatggaat ggctggctac
      421 tatggagtca gcctggccga ctgtaagtct cttcttgatg gctccagcta atcgctcaga
      481 acagggacag gggagcagtg agtgaccaga gaggaaactg cagcatacac tgtcatggtt
      541 gttccctgtg tgactgaggt catttgacag aaagcatcag tttttgtctt ttgaattaga
      601 gatgctgggg tgagggcttg gaagggacat tgtgtgccat ggtaaggcaa caggacctga
      661 ggggctcagg gactgtcact aaattttctg ttgcttcaga ctcaacatga gtaagggggc
      721 tgagaagaga tgggagagtc agccagagat cctggactga catggtattt gtcttcactt
      781 ctcttaccgt ggtgtgtaac ttttgagat
//
LOCUS       Lyz2-Cre|mutant          726 bp    dna     linear   UNK 
ACCESSION   replace
FEATURES             Location/Qualifiers
     Fragment        1..309
                     /name="left"
     Fragment        310..726
                     /name="change"
ORIGIN      
        1 tatttcacag cagcattgca gactagctaa aggcagaagg gagagactct ggagtgtctc
       61 aaatgtttaa gcattttaaa ataaatactg tatgcttatt cttgggctgc cagaatttct
      121 ctcatcacat aaatgaagaa ggaagatcaa gtgctgaagt ccatagatcg gtaggaactt
      181 cctgttttgc acacagctca aatgtaggaa accacaagct gttgggaaag gagggacttg
      241 gaggatgctt aaatagcagg catgctttct ctagtcagcc agcagctgac ccagcctcca
      301 gtcaccatgc ccaagaagaa gaggaaggtg tccaatttac tgaccgtaca ccaaaatttg
      361 cctgcattac cggtcgatgc aacgagtgat gaggttcgca agaacctgat ggacatgttc
      421 agggatcgcc aggcgttttc tgagcatacc tggaaaatgc ttctgtccgt ttgccggtcg
      481 tgggcggcat ggtgcaagtt gaataaccgg aaatggtttc ccgcagaacc tgaagatgtt
      541 cgcgattatc ttctatatct tcaggcgcgc ggtctggcag taaaaactat ccagcaacat
      601 ttgggccagc taaacatgct tcatcgtcgg tccgggctgc cacgaccaag tgacagcaat
      661 gctgtttcac tggttatgcg gcggatccga aaagaaaacg ttgatgccgg tgaacgtgca
      721 aaacag
//

=head2 Intermediate seq in FASTA format

>Lyz2-Cre|wild
TATTTCACAGCAGCATTGCAGACTAGCTAAAGGCAGAAGGGAGAGACTCTGGAGTGTCTC
AAATGTTTAAGCATTTTAAAATAAATACTGTATGCTTATTCTTGGGCTGCCAGAATTTCT
CTCATCACATAAATGAAGAAGGAAGATCAAGTGCTGAAGTCCATAGATCGGTAGGAACTT
CCTGTTTTGCACACAGCTCAAATGTAGGAAACCACAAGCTGTTGGGAAAGGAGGGACTTG
GAGGATGCTTAAATAGCAGGCATGCTTTCTCTAGTCAGCCAGCAGCTGACCCAGCCTCCA
GTCACCATGAAGACTCTCCTGACTCTGGGACTCCTCCTGCTTTCTGTCACTGCTCAGGCC
AAGGTCTATGAACGTTGTGAGTTTGCCAGAACTCTGAAAAGGAATGGAATGGCTGGCTAC
TATGGAGTCAGCCTGGCCGACTGTAAGTCTCTTCTTGATGGCTCCAGCTAATCGCTCAGA
ACAGGGACAGGGGAGCAGTGAGTGACCAGAGAGGAAACTGCAGCATACACTGTCATGGTT
GTTCCCTGTGTGACTGAGGTCATTTGACAGAAAGCATCAGTTTTTGTCTTTTGAATTAGA
GATGCTGGGGTGAGGGCTTGGAAGGGACATTGTGTGCCATGGTAAGGCAACAGGACCTGA
GGGGCTCAGGGACTGTCACTAAATTTTCTGTTGCTTCAGACTCAACATGAGTAAGGGGGC
TGAGAAGAGATGGGAGAGTCAGCCAGAGATCCTGGACTGACATGGTATTTGTCTTCACTT
CTCTTACCGTGGTGTGTAACTTTTGAGAT
>Lyz2-Cre|mutant
TATTTCACAGCAGCATTGCAGACTAGCTAAAGGCAGAAGGGAGAGACTCTGGAGTGTCTC
AAATGTTTAAGCATTTTAAAATAAATACTGTATGCTTATTCTTGGGCTGCCAGAATTTCT
CTCATCACATAAATGAAGAAGGAAGATCAAGTGCTGAAGTCCATAGATCGGTAGGAACTT
CCTGTTTTGCACACAGCTCAAATGTAGGAAACCACAAGCTGTTGGGAAAGGAGGGACTTG
GAGGATGCTTAAATAGCAGGCATGCTTTCTCTAGTCAGCCAGCAGCTGACCCAGCCTCCA
GTCACCATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTG
CCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCTGATGGACATGTTC
AGGGATCGCCAGGCGTTTTCTGAGCATACCTGGAAAATGCTTCTGTCCGTTTGCCGGTCG
TGGGCGGCATGGTGCAAGTTGAATAACCGGAAATGGTTTCCCGCAGAACCTGAAGATGTT
CGCGATTATCTTCTATATCTTCAGGCGCGCGGTCTGGCAGTAAAAACTATCCAGCAACAT
TTGGGCCAGCTAAACATGCTTCATCGTCGGTCCGGGCTGCCACGACCAAGTGACAGCAAT
GCTGTTTCACTGGTTATGCGGCGGATCCGAAAAGAAAACGTTGATGCCGGTGAACGTGCA
AAACAG

=head2 Intermediate SNP in FASTA format

>rs30925746,LFT,A/C
GCCCCTCAGGTCCTGTTGCCTTACCATGGCACACAATGTCCCTTCCAAGCCCTCACCCCAGCATCTCTAATTCAAAAGACAAAAACTGATGCTTTCTGTC
>rs30925746,RGT,A/C
AATGACCTCAGTCACACAGGGAACAACCATGACAGTGTATGCTGCAGTTTCCTCTCTGGTCACTCACTGCTCCCCTGTCCCTGTTCTGAGCGATTAGCTG

>rs30925748,LFT,C/T
CAATGTCCCTTCCAAGCCCTCACCCCAGCATCTCTAATTCAAAAGACAAAAACTGATGCTTTCTGTCAAATGACCTCAGTCACACAGGGAACAACCATGA
>rs30925748,RGT,C/T
AGTGTATGCTGCAGTTTCCTCTCTGGTCACTCACTGCTCCCCTGTCCCTGTTCTGAGCGATTAGCTGGAGCCATCAAGAAGAGACTTACAGTCGGCCAGG

>rs30925750,LFT,C/T
GCCCTCACCCCAGCATCTCTAATTCAAAAGACAAAAACTGATGCTTTCTGTCAAATGACCTCAGTCACACAGGGAACAACCATGACAGTGTATGCTGCAG
>rs30925750,RGT,C/T
TTCCTCTCTGGTCACTCACTGCTCCCCTGTCCCTGTTCTGAGCGATTAGCTGGAGCCATCAAGAAGAGACTTACAGTCGGCCAGGCTGACTCCATAGTAG

>rs30925752,LFT,C/G
ACCATGACAGTGTATGCTGCAGTTTCCTCTCTGGTCACTCACTGCTCCCCTGTCCCTGTTCTGAGCGATTAGCTGGAGCCATCAAGAAGAGACTTACAGT
>rs30925752,RGT,C/G
GGCCAGGCTGACTCCATAGTAGCCAGCCATTCCATTCCTTTTCAGAGTTCTGGCAAACTCACAACGTTCATAGACCTTGGCCTGAGCAGTGACAGAAAGC

>rs30926644,LFT,C/T
TTGGCCTGAGCAGTGACAGAAAGCAGGAGGAGTCCCAGAGTCAGGAGAGTCTTCATGGTGACTGGAGGCTGGGTCAGCTGCTGGCTGACTAGAGAAAGCA
>rs30926644,RGT,C/T
GCCTGCTATTTAAGCATCCTCCAAGTCCCTCCTTTCCCAACAGCTTGTGGTTTCCTACATTTGAGCTGTGTGCAAAACAGGAAGTTCCTACCGATCTATG

>rs30926646,LFT,G/T
AGGCTGGGTCAGCTGCTGGCTGACTAGAGAAAGCATGCCTGCTATTTAAGCATCCTCCAAGTCCCTCCTTTCCCAACAGCTTGTGGTTTCCTACATTTGA
>rs30926646,RGT,G/T
CTGTGTGCAAAACAGGAAGTTCCTACCGATCTATGGACTTCAGCACTTGATCTTCCTTCTTCATTTATGTGATGAGAGAAATTCTGGCAGCCCAAGAATA

=cut

=head1 Output Files
=head2 Annotated GenBank file

LOCUS       Lyz2-Cre|wild          809 bp    dna     linear   UNK 
DEFINITION  ORDER=1; TARGET=307,4
ACCESSION   replace
KEYWORDS    .
FEATURES             Location/Qualifiers
     Fragment        1..309
                     /name="left"
     Fragment        310..809
                     /name="change"
     SNP             197
                     /name="rs30926646"
                     /allele="G/T"
     SNP             262
                     /name="rs30926644"
                     /allele="C/T"
     SNP             439
                     /name="rs30925752"
                     /allele="C/G"
     SNP             517
                     /name="rs30925750"
                     /allele="C/T"
     SNP             532
                     /name="rs30925748"
                     /allele="C/T"
     SNP             565
                     /name="rs30925746"
                     /allele="A/C"
ORIGIN      
        1 tatttcacag cagcattgca gactagctaa aggcagaagg gagagactct ggagtgtctc
       61 aaatgtttaa gcattttaaa ataaatactg tatgcttatt cttgggctgc cagaatttct
      121 ctcatcacat aaatgaagaa ggaagatcaa gtgctgaagt ccatagatcg gtaggaactt
      181 cctgttttgc acacagctca aatgtaggaa accacaagct gttgggaaag gagggacttg
      241 gaggatgctt aaatagcagg catgctttct ctagtcagcc agcagctgac ccagcctcca
      301 gtcaccatga agactctcct gactctggga ctcctcctgc tttctgtcac tgctcaggcc
      361 aaggtctatg aacgttgtga gtttgccaga actctgaaaa ggaatggaat ggctggctac
      421 tatggagtca gcctggccga ctgtaagtct cttcttgatg gctccagcta atcgctcaga
      481 acagggacag gggagcagtg agtgaccaga gaggaaactg cagcatacac tgtcatggtt
      541 gttccctgtg tgactgaggt catttgacag aaagcatcag tttttgtctt ttgaattaga
      601 gatgctgggg tgagggcttg gaagggacat tgtgtgccat ggtaaggcaa caggacctga
      661 ggggctcagg gactgtcact aaattttctg ttgcttcaga ctcaacatga gtaagggggc
      721 tgagaagaga tgggagagtc agccagagat cctggactga catggtattt gtcttcactt
      781 ctcttaccgt ggtgtgtaac ttttgagat
//
LOCUS       Lyz2-Cre|mutant          726 bp    dna     linear   UNK 
DEFINITION  ORDER=0; TARGET=307,4
ACCESSION   replace
KEYWORDS    .
FEATURES             Location/Qualifiers
     Fragment        1..309
                     /name="left"
     Fragment        310..726
                     /name="change"
     SNP             197
                     /name="rs30926646"
                     /allele="G/T"
ORIGIN      
        1 tatttcacag cagcattgca gactagctaa aggcagaagg gagagactct ggagtgtctc
       61 aaatgtttaa gcattttaaa ataaatactg tatgcttatt cttgggctgc cagaatttct
      121 ctcatcacat aaatgaagaa ggaagatcaa gtgctgaagt ccatagatcg gtaggaactt
      181 cctgttttgc acacagctca aatgtaggaa accacaagct gttgggaaag gagggacttg
      241 gaggatgctt aaatagcagg catgctttct ctagtcagcc agcagctgac ccagcctcca
      301 gtcaccatgc ccaagaagaa gaggaaggtg tccaatttac tgaccgtaca ccaaaatttg
      361 cctgcattac cggtcgatgc aacgagtgat gaggttcgca agaacctgat ggacatgttc
      421 agggatcgcc aggcgttttc tgagcatacc tggaaaatgc ttctgtccgt ttgccggtcg
      481 tgggcggcat ggtgcaagtt gaataaccgg aaatggtttc ccgcagaacc tgaagatgtt
      541 cgcgattatc ttctatatct tcaggcgcgc ggtctggcag taaaaactat ccagcaacat
      601 ttgggccagc taaacatgct tcatcgtcgg tccgggctgc cacgaccaag tgacagcaat
      661 gctgtttcac tggttatgcg gcggatccga aaagaaaacg ttgatgccgg tgaacgtgca
      721 aaacag
//

=cut

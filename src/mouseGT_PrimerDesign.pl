#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  mouseGT_PrimerDesign.pl
#
#        USAGE:  ./mouseGT_PrimerDesign.pl -g gbfile -o prmout -d dir -t 1 -s 1
#
#  DESCRIPTION:  Mask SNPs and repeats in the GenBank sequence;
#                Design mouse genotyping primers
#
#      OPTIONS:  -g genbank file name: sequence file with fragments and SNP annotations 
#                -o primer design result file
#                -d temp directory
#                -t primer3 temperature option
#                   0 default (Breslauer KJ 1986 & Rychlik W 1990)
#                   1 recommend (SantaLucia JR 1998)
#                -s primer3 salt options
#                   0 default (Schildkraut C 1965)
#                   1 recommend (SantaLucia JR 1998)
#                   2 latest (Owczarzy R 2004)
# REQUIREMENTS:  Bioperl
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Mr. Hong Xu (hxu), hxu DOT hong AT gmail DOT com
#      COMPANY:  NIEHS
#      VERSION:  1.1
#      CREATED:  04/17/2010 10:51:40 AM
#     REVISION:  ---
#                1.1.0  (04/19/2012)
#                        @ get mutation type from seqobj->accession_number instead of seqobj->desc
#                        @ get primer design input from seqobi->desc
#                1.1.1  (06/23/2012)
#                        @ update primer3 version from 2.2.2 beta to 2.3.2 release (March 19, 2012)
#                1.1.2  (04/05/2013)
#                        @ write the dummy package - Primer3Pair to handle the problem of destroying 
#                          forward_primer or reverse_primer after accessing these methods in
#                          Bio::Tools::Primer3Redux::PrimerPair module (search 6/5/2012 hxu to see
#                          the changes) 
#                1.1.3  (05/25/2014)
#                        @ use Bio::Root->clone to deep clone Bio::Tools::Primer3Redux::Primer
#                1.1.4  (09/11/2014)
#                        @ break into two parts (flox & non-flox)
#                1.1.5  (10/03/2014)
#                        @ deep copy array of objects using Storable::dclone
#
#===============================================================================

use lib './modules'; # (hxu: To use the dummy Primer3Pair object)
use Primer3Pair; # (hxu: To use the dummy Primer3Pair object)
use strict;
use Getopt::Std;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::Run::Primer3Redux;
use Bio::Tools::Primer3Redux;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use File::Slurp;
use Storable qw(dclone);

############### Programming design #############
# initialize primer design region table (including "INCLUDED" and one of following:
#   "JUNCTION" - for 'deletion' in mutant, 'insertion' in wild
#   "TARGET" - for 'deletion' in wild, 'replace' in both wild and mutant
#   "POINT" - for 'point' in wild and mutant ('SEQUENCE_FORCE_LEFT_END' or 'SEQUENCE_FORCE_RIGHT_END')

my %prmins = ( 'INCLUDED' => 'SEQUENCE_INCLUDED_REGION',
               'JUNCTION' => 'SEQUENCE_OVERLAP_JUNCTION_LIST',
               'TARGET' => 'SEQUENCE_TARGET',
               'EXCLUDED' => 'SEQUENCE_EXCLUDED_REGION'
             );


# Tm calculation formula: 0, default (Breslauer KJ 1986 & Rychlik W 1990)
#                         1, recommend (SantaLucia JR 1998)
my $tm_formula = 0; 
# salt correction: 0, default (Schildkraut C 1965)
#                  1, recommend (SantaLucia JR 1998)
#                  2, latest (Owczarzy R 2004)
my $salt_correction = 0;

############### Get arguments #############
# get command line arguments
my %args;

getopt("gotsd", \%args);
if (!$args{o}) {
  die "Please specify primer output file name!\n";
}
else {
  open(FOT, ">$args{o}") or die "Cannot create primer output file - $args{o}";
}
if (!$args{g}) {
  die "Please specify GenBank file name!";
}
my $gbio = Bio::SeqIO->new(-file => "$args{g}",
                           -format => 'GenBank');
if ($args{t} == 1) {
  $tm_formula = 1;
}
if ( ($args{s} == 1) || ($args{s} == 2) ) {
  $salt_correction = $args{s};
}

my $dtop = $args{d};
my $ftmp = $dtop . '/' . md5_hex(rand().$$.time()); # set temp fasta file name
my $fai = $ftmp . '.fa';
my $famsk = $ftmp . '.mfa';
open(FLS, ">$ftmp.log2") or die "Cannot create log file - $ftmp.log2";

############### mask sequences #############
=head3 coding thoughts
   1. get seq features (fragment, SNP) by feature type when convert seqs
      from GenBank to FASTA
      1.1 Use fragment info for primer design parameters
      1.2 Use SNP info to mask SNP (soft)
   2. get mutation type when convert seqs from GenBank from FASTA
=cut

# convert seq from Genbank to FASTA; get features from seq object
# soft mask SNP
my $strnNm; # (hxu 07/12/2013: Add for mouse strain name)
my @sfeats = ();
my @sprms = ();
my $mutype;
my @nord = ();
my $flxsts;
while (my $seq = $gbio->next_seq) {
  # convert seq to FASTA format
  my $pid = $seq->display_id();
  $mutype = $seq->accession_number;
  my $prmdesc = $seq->desc;
  my @aprms = split('; ', $prmdesc);
  my $porder = shift @aprms;
  $porder =~ s/ORDER=//i;
  my ($sgen, $styp) = split('\|', $pid);
  $strnNm = $sgen unless $strnNm; # (hxu 07/12/2012: get strain name if strain name is unknown)
  my @sfs = $seq->get_SeqFeatures();
  my @frg = ();
  my @snp = ();
  # sort fragment & SNP features
  foreach my $feat (@sfs) {
    if ($feat->primary_tag eq 'Fragment') {
      push @frg, $feat;
    }
    elsif ($feat->primary_tag eq 'SNP') {
      push @snp, $feat->start;
    }
  }
  # soft mask SNP in seq
  my $seqstr = $seq->seq();
  $seqstr = strlc($seqstr, \@snp);
  $seq->seq($seqstr);
  
  # write seq, seq fragment features, and primer design parameters
  my $fai_n = $fai . '.' . $porder;
  my $faisq = Bio::SeqIO->new(-file => ">$fai_n",
                              -format => 'fasta');
  $faisq->write_seq($seq);
  $sfeats[$porder] = \@frg; # fragment features
  $sprms[$porder] = \@aprms;
  push @nord, $porder;
  $faisq->flush(); # flush file handle
  $faisq->DESTROY(); # close file handle

  # (hxu 09/04/2012: make deletion sequence file for flox mice)
  if ( ($mutype eq 'flox') && ($porder == 1) ) {
    my $fai_2 = $fai . '.' . 2;
    my $faisq2 = Bio::SeqIO->new(-file => ">$fai_2",
                              -format => 'fasta');

    # array of floxsties:
    #   1.before, 1.flox, 1.after, 2.before, 2.flox, 2.after
    my @newfrg = @{ dclone(\@frg) }; ##### (hxu 10/03/2012): deep copy array of objects
    my @newprm = @{ dclone(\@aprms) }; ##### (hxu 10/03/2012): deep copy array of objects
    (my $newsq, my $myflxfs, my $myprms, $flxsts) = getFloxsites($seq, \@newfrg, \@newprm);

    $faisq2->write_seq($newsq);
    $sfeats[2] = $myflxfs; # fragment features
    $sprms[2]  = $myprms;
    push @nord, 2;
    $faisq2->flush(); # flush file handle
    $faisq2->DESTROY(); # close file handle    
  }
}
$gbio->flush(); # flush file handle
$gbio->DESTROY(); # close file handle

# gmasker the sequence (soft masking)
foreach my $iOrd (@nord) {
  my $mfa_n = $famsk . '.' . $iOrd;
  my $fai_n = $fai . '.' . $iOrd;
  system("cat $fai_n | ./bin/gmasker -nbases 16 ./genome/glist/mouseRepeat.glist l both > $mfa_n");
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

##### break into two parts: one for flox, the other for non-flox
if ($mutype eq 'flox') {
  # flox mice
  # design universal energy-transfer-labeled primers
  my @uETLPrms_flx = labelFloxPrimer($mutype, $famsk, \@sprms, \@sfeats, 'flx');
  my @uETLPrms_mut = labelFloxPrimer($mutype, $famsk, \@sprms, \@sfeats, 'mut');

  # design melting curve primers
  my @meltPrms_flx = meltFloxPrimer($mutype, $famsk, \@sprms, \@sfeats, 'flx');
  my @meltPrms_mut = meltFloxPrimer($mutype, $famsk, \@sprms, \@sfeats, 'mut');

  # design size different primers
  my @sizePrms = sizeFloxPrimer($mutype, $famsk, \@sprms, \@sfeats);

  ############### Print primer #############
  # output primers

  my %outPrms = (
                  'Universal energy-transfer-labeled primers (Wild vs. Flox)' => \@uETLPrms_flx,
                  'Universal energy-transfer-labeled primers (Wild vs. Mut)'  => \@uETLPrms_mut,
                  'Melting curve primers (Wild vs. Flox)' => \@meltPrms_flx,
                  'Melting curve primers (Wild vs. Mut)'  => \@meltPrms_mut,                  
                  'Size different primers' => \@sizePrms,
                );
  prtFloxPrimers(\%outPrms, $strnNm);

  my $frm = $ftmp . '.' . '*';
  #system("rm $frm");
}
else {
  # non-flox mice
  # design universal energy-transfer-labeled primers
  my @uETLPrms = labelPrimer($mutype, $famsk, \@sprms, \@sfeats);

  # design melting curve primers
  my @meltPrms = meltPrimer($mutype, $famsk, \@sprms, \@sfeats);

  # design size different primers
  my @sizePrms = sizePrimer($mutype, $famsk, \@sprms, \@sfeats);

  ############### Print primer #############
  # output primers

  my %outPrms = (
                  'Universal energy-transfer-labeled primers' => \@uETLPrms,
                  'Melting curve primers' => \@meltPrms,
                  'Size different primers' => \@sizePrms,
                );
  prtPrimers(\%outPrms, $strnNm); # (hxu 07/12/2012: add $strnNm)

  my $frm = $ftmp . '.' . '*';
  # system("rm $frm");
}


################################################################################
################################################################################
#######                                                                  #######
#######                           sub routines                           #######
#######                                                                  #######
################################################################################
################################################################################

################################################################################
#######                                                                  #######
#######                         Main sub programs                        #######
#######                                                                  #######
################################################################################

# design energy-transfer-labeled primers
# design primers first, add tails to allele specific primers:
#   green for wild type; red for mutant
=head1 Universal Energy-Transfer-Labeled Primers
   Reference: Myakishev, M.V., Khripin, Y., Hu, S. & Hamer, D.H. 
              High-throughput SNP genotyping by allele-specific PCR with universal energy-transfer-labeled primers.
              Genome Res 11, 163-169 (2001).

        Fluoresein           Dabsyl (quencher)
        |(green)             |
        agcgatgcgttcgagcatcgctGAAGGTGACCAAGTTCATGCT
                              |
                              tail_1

        Sulforhodamine     Dabsyl (quencher)
        |(red)             |
        aggacgctgagatgcgtcctGAAGGTCGGAGTCAACGGATT
                            |
                            tail_2
=cut

sub labelPrimer {
  my ($mtype, $faf, $pp, $fts) = @_;

  my @results = ();
  
  my @flord = (0, 1); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

  if ($mtype eq 'point') {
    ##### design primer for 1st parameter set (left primer end or right primer end)
    # SEQUENCE_FORCE_LEFT_END
    # SEQUENCE_FORCE_RIGHT_END

    my $fp3m_0 = $faf . '.0'; 
 
    # set initial primer design parameters
    my %pmtr_0L = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula, 
                'PRIMER_PICK_ANYWAY' => 0,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '100-320',
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4, 
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);
    my %pmtr_0R = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula, 
                'PRIMER_PICK_ANYWAY' => 0,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[0]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      if ($k eq 'POINT') {
        my ($vend1, $vend2) = split('-', $v);
        $pmtr_0L{'SEQUENCE_FORCE_LEFT_END'} = $vend2 - 1;
        $pmtr_0R{'SEQUENCE_FORCE_RIGHT_END'} = $vend1 - 1;
      }
      else {         
        my $kv = $prmins{$k};
        $pmtr_0L{$kv} = $v;
        $pmtr_0R{$kv} = $v;
      }
    }

    # get allele info for seq_0
    my $allele_0;
    foreach my $ft (@{$fts->[0]}) {
      my @tagvalues = $ft->get_tag_values('name');
      if ($tagvalues[0] =~ m/allele$/) {
        my $seqio_0 = Bio::SeqIO->new(-file => $fp3m_0,
                                      -format => 'fasta');
        my $seq_0 = $seqio_0->next_seq;
        $allele_0 = $seq_0->subseq($ft->start, $ft->end);
        last;
      }
    }

    # design 1st primer pair set for both left & right end
    my @prm_0L = genPrimer($fp3m_0, \%pmtr_0L, $fts->[0]);
    my @prm_0R = genPrimer($fp3m_0, \%pmtr_0R, $fts->[0]);

    ##### design primer for 2nd parameter set
    # by changing the allele sub-sequence in 1st primer set
    my $fp3m_1 = $faf . '.1';

    # get allele info for seq_1
    my $allele_1;
    foreach my $ft (@{$fts->[1]}) {
      my @tagvalues = $ft->get_tag_values('name');
      if ($tagvalues[0] =~ m/allele$/) {
        my $seqio_1 = Bio::SeqIO->new(-file => $fp3m_1,
                                      -format => 'fasta');
        my $seq_1 = $seqio_1->next_seq;
        $allele_1 = $seq_1->subseq($ft->start, $ft->end);
        last;
      }
    }     

    # parse left end primer of 1st primer set to generate 2nd primer set
    my @prm_1L = ();
    foreach my $pp0L (@prm_0L) {
      my ($fp, $rp) = ($pp0L->forward_primer, $pp0L->reverse_primer);
      my $fp_1 = prmClone($fp);
      my $rp_1 = prmClone($rp);
      my $newseq = $fp_1->seq->seq;
      my $seqlen = length($newseq) - length($allele_1);
      $newseq = substr($newseq, 0, $seqlen) . $allele_1;
      $fp_1->seq->seq($newseq);
      my $ptm = primerTm($newseq);
      $fp_1->melting_temp($ptm);

      my $pp1L = new Primer3Pair; # (hxu 06/25/2012): dummy primer pair
      $pp1L->forward_primer($fp_1);
      $pp1L->reverse_primer($rp_1);
      push @prm_1L, $pp1L;  
    }
    
    # parse right end primer of 1st primer set to generate 2nd primer set
    my @prm_1R = ();
    foreach my $pp0R (@prm_0R) {
      my ($fp, $rp) = ($pp0R->forward_primer, $pp0R->reverse_primer);
      my $fp_1 = prmClone($fp);
      my $rp_1 = prmClone($rp);
      my $newseq = $rp_1->seq->seq;
      my $seqlen = length($allele_1) - length($newseq);
      $newseq = rvsCmp($allele_1) . substr($newseq, $seqlen);
      $rp_1->seq->seq($newseq);
      my $ptm = primerTm($newseq);
      $rp_1->melting_temp($ptm);

      my $pp1R = new Primer3Pair; # (hxu 06/25/2012): dummy primer pair
      $pp1R->forward_primer($fp_1);
      $pp1R->reverse_primer($rp_1);
      push @prm_1R, $pp1R;  
    }   

    # check primer compatibility
    my @chkResults_L = primerCompatible($faf, \@prm_0L, \@prm_1L, \@flord);
    push @results, @chkResults_L;
    my @chkResults_R = primerCompatible($faf, \@prm_0R, \@prm_1R, \@flord);
    push @results, @chkResults_R;
  }
  else {
    ##### design primer for 1st parameter set (junction)
    my $fp3m_0 = $faf . '.0'; 
 
    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[0]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_0{$kv} = $v;
    }

    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[0]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = 'fwd';
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = 'rvs';
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.1';
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '100-320',
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[1]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        $pmtr_1{$kv} = $v;
      }

      # add common primerd
      if ($cmnPT eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[1]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }
  }

  # add uETL tail
  # tail 1, 2 for universal energy-transfer-labeled primers
  my $utail1 = 'gaaggtgaccaagttcatgct'; # tail for green label, wild type
  my $utail2 = 'gaaggtcggagtcaacggatt'; # tail for red label, mutant
  foreach my $r (@results) {
    $r->[5] = $utail1 .'-'. $r->[5];
    $r->[9] = $utail2 .'-'. $r->[9];
  }  
  return @results;
}

# design melting curve primers
sub meltPrimer {
  my ($mtype, $faf, $pp, $fts) = @_;

  my @results = ();

  my @flord = (0, 1); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

  if ($mtype eq 'point') {
    ##### design primer for 1st parameter set (left primer end or right primer end)
    # SEQUENCE_FORCE_LEFT_END
    # SEQUENCE_FORCE_RIGHT_END

    my $fp3m_0 = $faf . '.0'; 
 
    # set initial primer design parameters
    my %pmtr_0L = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula, 
                'PRIMER_PICK_ANYWAY' => 0,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '80-280', 
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);
    my %pmtr_0R = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula, 
                'PRIMER_PICK_ANYWAY' => 0,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '80-280', 
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[0]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      if ($k eq 'POINT') {
        my ($vend1, $vend2) = split('-', $v);
        $pmtr_0L{'SEQUENCE_FORCE_LEFT_END'} = $vend2 - 1;
        $pmtr_0R{'SEQUENCE_FORCE_RIGHT_END'} = $vend1 - 1;
      }
      else {         
        my $kv = $prmins{$k};
        $pmtr_0L{$kv} = $v;
        $pmtr_0R{$kv} = $v;
      }
    }

    # get allele info for seq_0
    my $allele_0;
    foreach my $ft (@{$fts->[0]}) {
      my @tagvalues = $ft->get_tag_values('name');
      if ($tagvalues[0] =~ m/allele$/) {
        my $seqio_0 = Bio::SeqIO->new(-file => $fp3m_0,
                                      -format => 'fasta');
        my $seq_0 = $seqio_0->next_seq;
        $allele_0 = $seq_0->subseq($ft->start, $ft->end);
        last;
      }
    }

    # design 1st primer pair set for both left & right end
    my @prm_0L = genPrimer($fp3m_0, \%pmtr_0L, $fts->[0]);
    my @prm_0R = genPrimer($fp3m_0, \%pmtr_0R, $fts->[0]);

    ##### design primer for 2nd parameter set
    # by changing the allele sub-sequence in 1st primer set
    my $fp3m_1 = $faf . '.1';

    # get allele info for seq_1
    my $allele_1;
    foreach my $ft (@{$fts->[1]}) {
      my @tagvalues = $ft->get_tag_values('name');
      if ($tagvalues[0] =~ m/allele$/) {
        my $seqio_1 = Bio::SeqIO->new(-file => $fp3m_1,
                                      -format => 'fasta');
        my $seq_1 = $seqio_1->next_seq;
        $allele_1 = $seq_1->subseq($ft->start, $ft->end);
        last;
      }
    }     

    # parse left end primer of 1st primer set to generate 2nd primer set
    my @prm_1L = ();
    foreach my $pp0L (@prm_0L) {
      my ($fp, $rp) = ($pp0L->forward_primer, $pp0L->reverse_primer);
      my $fp_1 = prmClone($fp);
      my $rp_1 = prmClone($rp);
      my $newseq = $fp_1->seq->seq;
      my $seqlen = length($newseq) - length($allele_1);
      $newseq = substr($newseq, 0, $seqlen) . $allele_1;
      $fp_1->seq->seq($newseq);
      my $ptm = primerTm($newseq);
      $fp_1->melting_temp($ptm);

      my $pp1L = new Primer3Pair; # (hxu 06/25/2012): dummy primer pair
      $pp1L->forward_primer($fp_1);
      $pp1L->reverse_primer($rp_1);
      push @prm_1L, $pp1L;  
    }
    
    # parse right end primer of 1st primer set to generate 2nd primer set
    my @prm_1R = ();
    foreach my $pp0R (@prm_0R) {
      my ($fp, $rp) = ($pp0R->forward_primer, $pp0R->reverse_primer);
      my $fp_1 = prmClone($fp);
      my $rp_1 = prmClone($rp);
      my $newseq = $rp_1->seq->seq;
      my $seqlen = length($allele_1) - length($newseq);
      $newseq = rvsCmp($allele_1) . substr($newseq, $seqlen);
      $rp_1->seq->seq($newseq);
      my $ptm = primerTm($newseq);
      $rp_1->melting_temp($ptm);

      my $pp1R = new Primer3Pair;
      $pp1R->forward_primer($fp_1);
      $pp1R->reverse_primer($rp_1);
      push @prm_1R, $pp1R;  
    }   

    # check primer compatibility
    my @chkResults_L = primerCompatible($faf, \@prm_0L, \@prm_1L, \@flord);
    my @chkResults_R = primerCompatible($faf, \@prm_0R, \@prm_1R, \@flord);
    push @results, @chkResults_L;
    push @results, @chkResults_R;
  }
  else {
    ##### design primer for 1st parameter set (junction)
    my $fp3m_0 = $faf . '.0'; 
 
    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '80-230', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[0]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_0{$kv} = $v;
    }

    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[0]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.1';
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};
      
      # 2nd PCR product size limit
      my $tmprd_min = $cmnPT->[1] + 50;
      my $tmprd_max = $tmprd_min + 150;
      my $tmprd = $tmprd_min . '-' . $tmprd_max;
      
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => $tmprd,
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[1]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        $pmtr_1{$kv} = $v;
      }

      # add common primerd
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[1]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);

      push @results, @chkResults;
    }
  }

  ### add Tm shift tail
  #   if the two PCR amplicon Tms have bigger difference than 3°C, do not add tails
  #   otherwise add short tail to the low Tm primer, and add long tail to the high Tm primer
  # tail 1, 2 for Tm shift primers
  # my $tmtail1 = 'gataata'; # short tail for low Tm shift, used to be 'gcggac'
  my $tmtail1 = ' '; # (hxu 07/26/2012): do not add short tail per user's request
  my $tmtail2 = 'ggcagggcggcc'; # long tail for high Tm shift, used to be 'gcggcggcgggcagggcggcc'

  foreach my $r (@results) {
    my $Tmdif = $r->[8] - $r->[12];
    if (abs($Tmdif) < 3) {
      # check the common primer on forward or reverse strand
      my $cmp = $r->[3]; # common strand
      my $strand = 1;
      unless ($r->[13] =~ m/${cmp}/i) {
        # check strand on wild type product
        $strand = -1;
      }

      # add tail to primer and product
      if ($Tmdif < 0) {
        $r->[5] = $tmtail1 .'-'. $r->[5];
        $r->[9] = $tmtail2 .'-'. $r->[9];

        # if common primer on forward strand
        if ($strand > 0) {
          # wild product
          $r->[13] = $r->[13] . rvsCmp($tmtail1);
          # mutant product
          $r->[14] = $r->[14] . rvsCmp($tmtail2);
        }
        else {
          # wild product
          $r->[13] = $tmtail1 . $r->[13];
          # mutant product
          $r->[14] = $tmtail2 . $r->[14];;
        }
      }
      else {
        $r->[9] = $tmtail1 .'-'. $r->[9];
        $r->[5] = $tmtail2 .'-'. $r->[5];

        # if common primer on forward strand
        if ($strand > 0) {
          # wild product
          $r->[13] = $r->[13] . rvsCmp($tmtail2);
          # mutant product
          $r->[14] = $r->[14] . rvsCmp($tmtail1);
        }
        else {
          # wild product
          $r->[13] = $tmtail2 . $r->[13];
          # mutant product
          $r->[14] = $tmtail1 . $r->[14];;
        }
      }
      $r->[8]  = meltSim($r->[13]);
      $r->[12]  = meltSim($r->[14]);
    }
  }  
  return @results;
}

# design amplicon size different primers
sub sizePrimer {
  my ($mtype, $faf, $pp, $fts) = @_;

  my @results = ();

  my @flord = (0, 1); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

  if ($mtype eq 'point') {
    # size different primers not good for point mutation
    return @results;
  }
  else {
    ##### design primer for 1st parameter set (junction)
    my $fp3m_0 = $faf . '.0'; 
 
    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '100-250', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[0]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_0{$kv} = $v;
    }

    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[0]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.1';
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};
      
      # 2nd PCR product size limit
      my $tmprd_min = $cmnPT->[1] + 50;
      my $tmprd_max = $tmprd_min + 150;
      my $tmprd = $tmprd_min . '-' . $tmprd_max;
    
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => $tmprd, 
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[1]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        $pmtr_1{$kv} = $v;
      }


      # add common primerd
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[1]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }
  }

  return @results;
}

################################################################################
##### Flox mice primer design --> BEGIN
################################################################################

sub labelFloxPrimer {
  my ($mtype, $faf, $pp, $fts, $myvs) = @_;

  my @results = ();

  # use the global variable of $flxsts for floxed regions:
  #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after
  my $flen1 = $flxsts->[2]->[1] - $flxsts->[0]->[0];
  my $flen2 = $flxsts->[5]->[1] - $flxsts->[3]->[0];
  
  if ($myvs eq 'flx') {
  
    my @flord = (1, 0); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

    # wild vs. flox
    # design flox primer first on bigger floxsite
    # find the floxed site with the longer insertion
    # default $flen1 <= $flen2
    my $tgtv = ($flxsts->[5]->[1] - 3) . ',' . 4;
    if ($flen1 > $flen2) {
      $tgtv = ($flxsts->[0]->[0] - 3) . ',' . 4;
    }

    my $fp3m_0 = $faf . '.' . $flord[0]; 

    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                  'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                  'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                  'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                  'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[$flord[0]]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_0{$kv} = $v;
    }

    # reset the target 
    $pmtr_0{'SEQUENCE_TARGET'} = $tgtv;
      
    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[$flord[0]]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.' . $flord[1];
  
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};
    
      # design parameter of 2nd primer pair 
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                    'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                    'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 
                    'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                    'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[$flord[1]]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        $pmtr_1{$kv} = $v;
      }

      # add common primer
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[$flord[1]]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }
  }
  elsif ($myvs eq 'mut') {
    # wild vs. mutant

    my @flord = (0, 2); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

    # design wild primer first on junctions
    # use the global variable of $flxsts for floxed regions:
    #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after

    my $fp3m_0 = $faf . '.' . $flord[0]; 

    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                  'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                  'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                  'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                  'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[$flord[0]]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      next if ($kv eq 'SEQUENCE_TARGET'); # ignore target region limit
      $pmtr_0{$kv} = $v;
    }

    # design wild primer first on junctions
    # use the global variable of $flxsts for floxed regions:
    #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after
    
    $pmtr_0{'SEQUENCE_OVERLAP_JUNCTION_LIST'} = ($flxsts->[0]->[0] - 1) . ' ' . 
                                                ($flxsts->[0]->[0] + ($flxsts->[3]->[0] - $flxsts->[2]->[1]));
    $pmtr_0{'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'} = 3;
    $pmtr_0{'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'} = 10;
      
    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[$flord[0]]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.' . $flord[1];
  
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};
    
      # design parameter of 2nd primer pair 
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                    'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                    'PRIMER_PRODUCT_SIZE_RANGE' => '100-320', 
                    'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                    'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[$flord[1]]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        next if ($kv eq 'SEQUENCE_TARGET'); # ignore target region limit
        $pmtr_1{$kv} = $v;
      }

      # add common primer
      my $regionSta = $flxsts->[0]->[0] - 1;
      my $regionLen = ($flxsts->[0]->[1] - $flxsts->[0]->[0]) + 1 +
                      ($flxsts->[5]->[1] - $flxsts->[4]->[0]) + 1;
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
        # limit right primer picking region in mutant insertion
        $pmtr_1{'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'} = ',,' . $regionSta . ',' . $regionLen;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
        # limit right primer picking region in mutant insertion
        $pmtr_1{'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'} = $regionSta . ',' . $regionLen . ',,';
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[$flord[1]]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }    
  }

  # add uETL tail
  # tail 1, 2 for universal energy-transfer-labeled primers
  my $utail1 = 'gaaggtgaccaagttcatgct'; # tail for green label, wild type
  my $utail2 = 'gaaggtcggagtcaacggatt'; # tail for red label, mutant
  foreach my $r (@results) {
    $r->[5] = $utail1 .'-'. $r->[5];
    $r->[9] = $utail2 .'-'. $r->[9];
  }
  
  return @results;
}

# design melting curve primers
sub meltFloxPrimer {
  my ($mtype, $faf, $pp, $fts, $myvs) = @_;

  my @results = ();

  # use the global variable of $flxsts for floxed regions:
  #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after
  my $flen1 = $flxsts->[2]->[1] - $flxsts->[0]->[0];
  my $flen2 = $flxsts->[5]->[1] - $flxsts->[3]->[0];

  if ($myvs eq 'flx') {
  
    my @flord = (1, 0); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

    # wild vs. flox
    # design flox primer first on bigger floxsite
    # find the floxed site with the longer insertion
    # default $flen1 <= $flen2
    my $tgtv = ($flxsts->[5]->[1] - 3) . ',' . 4;
    if ($flen1 > $flen2) {
      $tgtv = ($flxsts->[0]->[0] - 3) . ',' . 4;
    }

    my $fp3m_0 = $faf . '.' . $flord[0]; 

    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                  'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                  'PRIMER_PRODUCT_SIZE_RANGE' => '80-230', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                  'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                  'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[$flord[0]]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_0{$kv} = $v;
    }

    # reset the target 
    $pmtr_0{'SEQUENCE_TARGET'} = $tgtv;
      
    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[$flord[0]]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.' . $flord[1];
  
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};

      # 2nd PCR product size limit
      my $tmprd_min = $cmnPT->[1] + 50;
      my $tmprd_max = $tmprd_min + 250;
      my $tmprd = $tmprd_min . '-' . $tmprd_max;
    
      # design parameter of 2nd primer pair 
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                    'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                    'PRIMER_PRODUCT_SIZE_RANGE' => $tmprd, 
                    'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                    'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[$flord[1]]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        $pmtr_1{$kv} = $v;
      }

      # add common primer
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[$flord[1]]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }
  }
  elsif ($myvs eq 'mut') {
    # wild vs. mutant

    my @flord = (0, 2); # (hxu 09/12/2012: add file order for compatible with flox primers - labeled & melt)

    # design wild primer first on junctions
    # use the global variable of $flxsts for floxed regions:
    #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after

    my $fp3m_0 = $faf . '.' . $flord[0]; 

    # set initial primer design parameters
    my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                  'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                  'PRIMER_PRODUCT_SIZE_RANGE' => '80-230', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                  'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                  'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[$flord[0]]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      next if ($kv eq 'SEQUENCE_TARGET'); # ignore target region limit
      $pmtr_0{$kv} = $v;
    }

    # design wild primer first on junctions
    # use the global variable of $flxsts for floxed regions:
    #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after
    
    $pmtr_0{'SEQUENCE_OVERLAP_JUNCTION_LIST'} = ($flxsts->[0]->[0] - 1) . ' ' . 
                                                ($flxsts->[0]->[0] + ($flxsts->[3]->[0] - $flxsts->[2]->[1]));
    $pmtr_0{'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'} = 3;
    $pmtr_0{'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION'} = 10;
      
    # design 1st primer pair set
    my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[$flord[0]]);

    # parse primer result to get common primer
    my %cmnPrms = ();
    foreach my $pair0 (@prm_0) {
      my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
      my $prdsize = $rp->end - $fp->start + 1; # PCR product size
      # check whether forward primer is common primer
      if ( $fp->oligo_type eq 'common' ) {
        my $sk = $fp->seq->seq;
        $cmnPrms{$sk} = ['fwd', $prdsize];
      }
      # check whether reverse primer is common primer
      elsif ( $rp->oligo_type eq 'common' ) {
        my $sk = $rp->seq->seq;
        $cmnPrms{$sk} = ['rvs', $prdsize];
      }
    }
  
    ##### design primer for 2nd parameter set
    my $fp3m_1 = $faf . '.' . $flord[1];
  
    foreach my $f1stseq (keys %cmnPrms) {
      my $cmnPT = $cmnPrms{$f1stseq};

      # 2nd PCR product size limit
      my $tmprd_min = $cmnPT->[1] + 50;
      my $tmprd_max = $tmprd_min + 250;
      my $tmprd = $tmprd_min . '-' . $tmprd_max;
    
      # design parameter of 2nd primer pair 
      my %pmtr_1 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                    'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                    'PRIMER_PRODUCT_SIZE_RANGE' => $tmprd, 
                    'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                    'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

      # parse specific primer design parameters
      foreach my $dp (@{$pp->[$flord[1]]}) {
        $dp =~ m/(.+)\=(.+)/;
        my $k = $1;
        my $v = $2;
        my $kv = $prmins{$k};
        next if ($kv eq 'SEQUENCE_TARGET'); # ignore target region limit
        $pmtr_1{$kv} = $v;
      }

      # add common primer
      my $regionSta = $flxsts->[0]->[0] - 1;
      my $regionLen = ($flxsts->[0]->[1] - $flxsts->[0]->[0]) + 1 +
                      ($flxsts->[5]->[1] - $flxsts->[4]->[0]) + 1;
      if ($cmnPT->[0] eq 'fwd') {
        $pmtr_1{'SEQUENCE_PRIMER'} = $f1stseq;
        # limit right primer picking region in mutant insertion
        $pmtr_1{'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'} = ',,' . $regionSta . ',' . $regionLen;
      }
      else {
        $pmtr_1{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
        # limit right primer picking region in mutant insertion
        $pmtr_1{'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'} = $regionSta . ',' . $regionLen . ',,';
      }

      # call general primer design sub module
      my @prm_1 = genPrimer($fp3m_1, \%pmtr_1, $fts->[$flord[1]]);

      # check primer compatibility
      my @chkResults = primerCompatible($faf, \@prm_0, \@prm_1, \@flord);
      push @results, @chkResults;
    }    
  }

  ### add Tm shift tail
  #   if the two PCR amplicon Tms have bigger difference than 3°C, do not add tails
  #   otherwise add short tail to the low Tm primer, and add long tail to the high Tm primer
  # tail 1, 2 for Tm shift primers
  # my $tmtail1 = 'gataata'; # short tail for low Tm shift, used to be 'gcggac'
  my $tmtail1 = ' '; # (hxu 07/26/2012): do not add short tail per user's request
  my $tmtail2 = 'ggcagggcggcc'; # long tail for high Tm shift, used to be 'gcggcggcgggcagggcggcc'

  foreach my $r (@results) {
    my $Tmdif = $r->[8] - $r->[12];
    if (abs($Tmdif) < 3) {
      # check the common primer on forward or reverse strand
      my $cmp = $r->[3]; # common strand
      my $strand = 1;
      unless ($r->[13] =~ m/${cmp}/i) {
        # check strand on wild type product
        $strand = -1;
      }

      # add tail to primer and product
      if ($Tmdif < 0) {
        $r->[5] = $tmtail1 .'-'. $r->[5];
        $r->[9] = $tmtail2 .'-'. $r->[9];

        # if common primer on forward strand
        if ($strand > 0) {
          # wild product
          $r->[13] = $r->[13] . rvsCmp($tmtail1);
          # mutant product
          $r->[14] = $r->[14] . rvsCmp($tmtail2);
        }
        else {
          # wild product
          $r->[13] = $tmtail1 . $r->[13];
          # mutant product
          $r->[14] = $tmtail2 . $r->[14];;
        }
      }
      else {
        $r->[9] = $tmtail1 .'-'. $r->[9];
        $r->[5] = $tmtail2 .'-'. $r->[5];

        # if common primer on forward strand
        if ($strand > 0) {
          # wild product
          $r->[13] = $r->[13] . rvsCmp($tmtail2);
          # mutant product
          $r->[14] = $r->[14] . rvsCmp($tmtail1);
        }
        else {
          # wild product
          $r->[13] = $tmtail2 . $r->[13];
          # mutant product
          $r->[14] = $tmtail1 . $r->[14];;
        }
      }
      $r->[8]  = meltSim($r->[13]);
      $r->[12]  = meltSim($r->[14]);
    }
  }
  
  return @results;
}

# design amplicon size different primers
sub sizeFloxPrimer {
  my ($mtype, $faf, $pp, $fts) = @_;
  my @results = ();

  # (hxu 08/07/12: added for flox-cre condition knockout mice)

  # There are three seuquences for consideration: wild, floxed, deleted
  # deleted can be generated from floxed (see sub getFloxsites())

  # design primer for wild type mice sequence first
  # use the global variable of $flxsts for floxed regions:
  #   floxsite.1.before, floxsite.1.flox, floxsite.1.after, floxsite.2.before, floxsite.2.flox, floxsite.2.after
  # find the floxed site with the shorter insertion
  my $flen1 = $flxsts->[2]->[1] - $flxsts->[0]->[0];
  my $flen2 = $flxsts->[5]->[1] - $flxsts->[3]->[0];
  
  # default $flen1 <= $flen2
  my $tgtv = ($flxsts->[0]->[0] - 10) . ',' . 20;
  if ($flen1 > $flen2) {
    $tgtv = (($flxsts->[0]->[0] - 1)
          + ($flxsts->[3]->[0] - $flxsts->[2]->[1]) - 10) . ',' . 20;
  }

  my $fp3m_0 = $faf . '.0'; 

  # set initial primer design parameters
  my %pmtr_0 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                'PRIMER_PRODUCT_SIZE_RANGE' => '60-220', 'PRIMER_MIN_THREE_PRIME_DISTANCE' => 0,
                'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

  # parse specific primer design parameters
  foreach my $dp (@{$pp->[0]}) {
    $dp =~ m/(.+)\=(.+)/;
    my $k = $1;
    my $v = $2;
    my $kv = $prmins{$k};
    $pmtr_0{$kv} = $v;
  }

  # reset the target for wild primer
  $pmtr_0{'SEQUENCE_TARGET'} = $tgtv;
      
  # design 1st primer pair set
  my @prm_0 = genPrimer($fp3m_0, \%pmtr_0, $fts->[0]);

  # parse primer result to get common primer
  my %cmnPrms = ();
  foreach my $pair0 (@prm_0) {
    my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
    my $prdsize = $rp->end - $fp->start + 1; # PCR product size
    # check whether forward primer is common primer
    if ( $fp->oligo_type eq 'common' ) {
      my $sk = $fp->seq->seq;
      $cmnPrms{$sk} = ['fwd', $prdsize];
    }
    # check whether reverse primer is common primer
    elsif ( $rp->oligo_type eq 'common' ) {
      my $sk = $rp->seq->seq;
      $cmnPrms{$sk} = ['rvs', $prdsize];
    }
  }
  
  ##### design primer for 2nd parameter set
  my $fp3m_2 = $faf . '.2';
  
  foreach my $f1stseq (keys %cmnPrms) {
    my $cmnPT = $cmnPrms{$f1stseq};

    # 2nd PCR product size limit
    my $tmprd_min = $cmnPT->[1] + $flen1 + 40;
    my $tmprd_max = $tmprd_min + 150;
    my $tmprd = $tmprd_min . '-' . $tmprd_max;
    
    # design parameter of 2nd primer pair 
    my %pmtr_2 = ('PRIMER_NUM_RETURN' => 10, 'PRIMER_LOWERCASE_MASKING' => 1, 'PRIMER_TM_FORMULA' => $tm_formula,
                  'PRIMER_SALT_CORRECTIONS' => $salt_correction, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH' => './config/primer3_config/',
                  'PRIMER_PRODUCT_SIZE_RANGE' => $tmprd, 
                  'PRIMER_OPT_SIZE' => 22, 'PRIMER_MAX_SIZE' => 30, 'PRIMER_MAX_END_GC' => 3, 'PRIMER_MAX_POLY_X' => 4,
                  'PRIMER_OPT_TM' => 64, 'PRIMER_MAX_TM' => 67, 'PRIMER_MIN_TM' => 60);

    # parse specific primer design parameters
    foreach my $dp (@{$pp->[2]}) {
      $dp =~ m/(.+)\=(.+)/;
      my $k = $1;
      my $v = $2;
      my $kv = $prmins{$k};
      $pmtr_2{$kv} = $v;
    }

    # add common primer
    if ($cmnPT->[0] eq 'fwd') {
      $pmtr_2{'SEQUENCE_PRIMER'} = $f1stseq;
    }
    else {
      $pmtr_2{'SEQUENCE_PRIMER_REVCOMP'} = $f1stseq;
    }

    # call general primer design sub module
    my @prm_2 = genPrimer($fp3m_2, \%pmtr_2, $fts->[2]);

    # check primer compatibility
    my @chkResults = primerSizeFloxCompatible($faf, \@prm_0, \@prm_2);
    push @results, @chkResults;
  }

  return @results;
}

################################################################################
##### Flox mice primer design <-- END
################################################################################


################################################################################
#######                                                                  #######
#######                    Supplementary sub programs                    #######
#######                                                                  #######
################################################################################

### change string chars to lower case (soft masking)
sub strlc {
  my ($seq, $loc) = @_;
  foreach my $lcn (@$loc) {
    my $char = substr($seq, $lcn - 1, 1);
    substr($seq, $lcn - 1, 1) = lc($char);
  }
  return $seq;
}

### melting Tm simulation program for PCR amplicon
# run meltSim and parse the result
sub meltSim {
  my $seq = shift;
  my $ftsq = $dtop . '/' . 'temp.meltsim.seq';
  my $fttm = $dtop . '/' . 'temp.meltsim.tm';
  open(FTSQ, ">$ftsq") or die "Cannot create seq file - $ftsq";
  print FTSQ $seq, "\n";
  close(FTSQ);
  system("cat $ftsq | ./bin/meltsim -g $fttm -o /dev/null"); # get Tm prediction
  my $stm = read_file($fttm);
  $stm =~ /Melting Temperature \(Tm\)  : (.+) /; 
  return $1;
}


### primer melting Tm
# run oligotm to get short oligo Tm xuh2
sub primerTm {
  my $oligo = shift;
  my $tm = `oligotm -tp 1 -sc 1 $oligo`;
  chomp($tm);
  return $tm;
}


##### design general primers
sub genPrimer {
  my ($fain, $prmds, $ft) = @_;

  # get strain info: wild or mutant from seq primary ID
  my $fseq = Bio::SeqIO->new(-file => "$fain",
                             -format => 'fasta');
  my $seqobj = $fseq->next_seq;
  my $seqid = $seqobj->display_id;
  my ($sgene, $strain) = split('\|', $seqid);
    
  my $pr3 = Bio::Tools::Run::Primer3Redux->new(-outfile => $dtop . '/' . 'gp3.temp.out',
                                               -path => "./bin/primer3_core");
  unless ($pr3->executable) {
    print FLS "primer3 can not be found. Is it installed?\n";
    exit(-1);
  }

  # set initial primer design parameters
  $pr3->set_parameters( %$prmds );

  my $seqio = Bio::SeqIO->new(-file => $fain,
                              -format => 'fasta');
  my $seq = $seqio->next_seq;
    
  # run primer3
  $pr3->run($seq);
  print FLS "Run primer3 ......\n";

  # parse primer to get unique common primers in the 1st primer set
  my $p3 = Bio::Tools::Primer3Redux->new( -file => $dtop . '/' . 'gp3.temp.out' );
  my @prmpairs = ();
  while (my $p = $p3->next_result) {
    while (my $pair = $p->next_primer_pair) {
      my ($fp, $rp) = ($pair->forward_primer, $pair->reverse_primer);
      $fp->oligo_type($strain);
      $rp->oligo_type($strain);

      # check common primer
      isCommonPrimer($fp, $ft);
      isCommonPrimer($rp, $ft);

      # (6/5/2012 hxu: to put the primer back to a dummy Primer3Pair object)
      my $p3p = new Primer3Pair;
      $p3p->forward_primer($fp);
      $p3p->reverse_primer($rp);
      # (6/5/2012 hxu: end).

      push @prmpairs, $p3p;
    }
  }
  print FLS "Parse primer3 ......\n";
  return @prmpairs;
}  
  

##### check whether primer is in common region
sub isCommonPrimer {
  my ($prm, $feats) = @_;
  foreach my $ft (@$feats) {
    my @tagvalues = $ft->get_tag_values('name');
    if ($tagvalues[0] =~ m/(left|right)$/i) {
      if ( ($prm->start > $ft->start) && ($prm->end < $ft->end) ) {
        $prm->oligo_type('common');
        return 1;
      }
    }
  }
  return 0;
}

##### check primer compatibilty
# each primer trio will be grouped by unique common primer
sub primerCompatible {
  my ($fa, $prm0, $prm1, $fext) = @_;

  # get seq objects
  my $fain0 = $fa . '.' . $fext->[0];
  my $seqio0 = Bio::SeqIO->new(-file => $fain0,
                               -format => 'fasta');
  my $seqobj0 = $seqio0->next_seq;
  my $fain1 = $fa . '.' . $fext->[1];
  my $seqio1 = Bio::SeqIO->new(-file => $fain1,
                               -format => 'fasta');
  my $seqobj1 = $seqio1->next_seq;

  # parse 1st primer set to get unique common primers
  my %cmnSet = ();
  foreach my $pair0 (@$prm0) {
    my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
    # check forward primer first
    if ($fp->oligo_type eq 'common') {
      my $seqstr = $fp->seq->seq;
      unless ( exists $cmnSet{$seqstr} ) {
        my @pps = ([], []);
        $cmnSet{$seqstr} = \@pps;
      }
      push @{$cmnSet{$seqstr}->[0]}, $pair0;
    }
    # then check reverse primer
    elsif ($rp->oligo_type eq 'common') {
      my $seqstr = $rp->seq->seq;
      unless ( exists $cmnSet{$seqstr} ) {
        my @pps = ([], []);
        $cmnSet{$seqstr} = \@pps;
      }
      push @{$cmnSet{$seqstr}->[0]}, $pair0;
    }
  }

  # parse 2nd primer set to group       
  foreach my $pair1 (@$prm1) {
    my ($fp, $rp) = ($pair1->forward_primer, $pair1->reverse_primer);
    my $fpseq = $fp->seq->seq;
    my $rpseq = $rp->seq->seq;

    # check forward primer first
    if (exists $cmnSet{$fpseq}) {
      push @{$cmnSet{$fpseq}->[1]}, $pair1;
    }
    elsif (exists $cmnSet{$rpseq}) {
      push @{$cmnSet{$rpseq}->[1]}, $pair1;
    }
  }

  # make combinatory primer sets and check multiplex scores
  my @prmResults = ();
  foreach my $sk (keys %cmnSet) {
    foreach my $ppr0 (@{$cmnSet{$sk}->[0]}) {
      foreach my $ppr1 (@{$cmnSet{$sk}->[1]}) {
        my @mscore = multiplxScore($ppr0, $seqobj0, $ppr1, $seqobj1);
        my $strn0 = shift(@mscore);
        my $strn1 = shift(@mscore);

        # prepare primer info in format:
        #   PrimPrimEnd2Score,PrimPrimEnd1Score,PrimPrimAnyScore
        #   CommPrimer,CommTm, WildPrimer,WildTm,WildProdSize,WildProdTm, MutPrimer,MutTm,MutProdSize,MutProdTm
        #   WildProdSeq,MutantProdSeq
        my @prmSet = convPrimerSet($sk, $seqobj0, $ppr0, $strn0, $seqobj1, $ppr1, $strn1);
        push(@mscore, @prmSet);
        push @prmResults, \@mscore;
      }
    }
  }

  return @prmResults;
}  

##### calculate MultiPLX scores
# higher is better (delta G is higher)
sub multiplxScore {
  my ($prm1, $seq1, $prm2, $seq2) = @_;
  my $fmp = $dtop . '/' . 'multiPLX.in';
  my $fms = $dtop . '/' . 'multiPLX.score';

  # prepare input file for MultiPLX
  open(FMP, ">$fmp") or die "Cannot open MultiPLX input file - $fmp";
  my $subseq1 = $seq1->subseq( $prm1->forward_primer->start, $prm1->reverse_primer->end );
  my $seqid1 = $seq1->display_id;
  my ($sgene1, $strain1) = split('\|', $seqid1);
  my $subseq2 = $seq2->subseq( $prm2->forward_primer->start, $prm2->reverse_primer->end );
  my $seqid2 = $seq2->display_id;
  my ($sgene2, $strain2) = split('\|', $seqid2);
  print FMP join("\t", $strain1, $prm1->forward_primer->seq->seq, $prm1->reverse_primer->seq->seq, $subseq1), "\n";
  print FMP join("\t", $strain2, $prm2->forward_primer->seq->seq, $prm2->reverse_primer->seq->seq, $subseq2), "\n";
  close(FMP);

  # run MultiPLX
  system("./bin/cmultiplx -thermodynamics ./config/thermodynamics.txt -primers $fmp -calcscores 123 -savescores $fms");

  # parse score file
  open(FMS, "<$fms") or die "Cannot open MultiPLX score file - $fms";
  <FMS>; # skip file header
  my $ln = <FMS>;
  chomp($ln);
  $ln =~ s/\x0D//; # remove Ctrl-M from end
  my @scores = split(/\t/, $ln);
  return @scores;  
}

# primer input set format: seqObj_0, primerObj_0, strain_0, seqObj_1, primerObj_1, strain_1
# primer output set format: CommPrimer,CommTm, WildPrimer,WildTm,WildProdSize,WildProdTm, MutPrimer,MutTm,MutProdSize,MutProdTm
sub convPrimerSet {
  my ($cmnSeq, $sqObj0, $pr0, $str0, $sqObj1, $pr1, $str1) = @_;
  my @result = ($cmnSeq);
  my ($fp0, $rp0) = ($pr0->forward_primer, $pr0->reverse_primer);
  my ($fp1, $rp1) = ($pr1->forward_primer, $pr1->reverse_primer);

  # check common primer
  if ($fp0->seq->seq eq $cmnSeq) {
    $result[1] = $fp0->melting_temp;
    if ($str0 eq 'wild') {
      # primer pair 0 gives wild fragment
      $result[2] = $rp0->seq->seq;
      $result[3] = $rp0->melting_temp;
      my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
      $result[4] = length($seqStr0);
      $result[5] = meltSim($seqStr0);
      # primer pair 1 gives mutant fragment
      $result[6] = $rp1->seq->seq;
      $result[7] = $rp1->melting_temp;
      my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
      $result[8] = length($seqStr1);
      $result[9] = meltSim($seqStr1);
      $result[10] = $seqStr0;
      $result[11] = $seqStr1;
    }
    else {
      # primer pair 0 gives mutant fragment
      $result[6] = $rp0->seq->seq;
      $result[7] = $rp0->melting_temp;
      my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
      $result[8] = length($seqStr0);
      $result[9] = meltSim($seqStr0);
      # primer pair 1 gives wild fragment
      $result[2] = $rp1->seq->seq;
      $result[3] = $rp1->melting_temp;
      my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
      $result[4] = length($seqStr1);
      $result[5] = meltSim($seqStr1);
      $result[10] = $seqStr1;
      $result[11] = $seqStr0;
    }
  }
  else {
    $result[1] = $rp0->melting_temp;
    if ($str0 eq 'wild') {
      # primer pair 0 gives wild fragment
      $result[2] = $fp0->seq->seq;
      $result[3] = $fp0->melting_temp;
      my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
      $result[4] = length($seqStr0);
      $result[5] = meltSim($seqStr0);
      # primer pair 1 gives mutant fragment
      $result[6] = $fp1->seq->seq;
      $result[7] = $fp1->melting_temp;
      my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
      $result[8] = length($seqStr1);
      $result[9] = meltSim($seqStr1);
      $result[10] = $seqStr0;
      $result[11] = $seqStr1;
    }
    else {
      # primer pair 0 gives mutant fragment
      $result[6] = $fp0->seq->seq;
      $result[7] = $fp0->melting_temp;
      my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
      $result[8] = length($seqStr0);
      $result[9] = meltSim($seqStr0);
      # primer pair 1 gives wild fragment
      $result[2] = $fp1->seq->seq;
      $result[3] = $fp1->melting_temp;
      my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
      $result[4] = length($seqStr1);
      $result[5] = meltSim($seqStr1);
      $result[10] = $seqStr1;
      $result[11] = $seqStr0;
    }
  }
  return @result;
}

 
# print out designed primers
sub prtPrimers {
  my $hPrms = shift;
  my $mmNm = shift; # (hxu 07/12/2012: add $mmNm for mouse strain name)

  print FOT 'Mouse Genotyping Primers for Mouse Strain ~ ' . $mmNm, "\n\n";
  
  foreach my $k (sort keys %$hPrms) {
    my $prms = $hPrms->{$k};
    print FOT 'Primer set - ' . $k . ':', "\n";
    print FOT join("\t", 'PrimPrimEnd2','PrimPrimEnd1','PrimPrimAny',
                         'CommPrimer','CommTm', 'WildPrimer','WildTm','WildProdSize','WildProdTm', 
                         'MutPrimer','MutTm','MutProdSize','MutProdTm', 'WildProdSeq','MutantProdSeq'), "\n";
    my @sorted = sort { $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] } @$prms;
    foreach my $pt (@sorted) {
      print FOT join("\t", @$pt), "\n";
    }
    print FOT "\n\n";
  }      
}

# primer object clone
sub prmClone {
  my $pin = shift;
  my $newprm = $pin->clone; # (hxu 06/25/2012): update with Bio::Root->clone
  return $newprm;
} 

# reverse complement
sub rvsCmp {
  my $sin = shift;
  return $sin unless ($sin =~ m/[ACGT]+/i);
  my $sqin = Bio::Seq->new( -seq => $sin );
  my $sqot = $sqin->revcom();  
  return $sqot->seq;
} 


################################################################################
##### Flox mice primer design --> BEGIN
################################################################################

##### isPcr search primers
sub primerFloxIsPcr {
  my ($fa, $prm) = @_;

  # prepare the isPcr search
  my $outfmt = '-out=bed';
  my $ptmp = $ftmp . 'prm.qry';
  my $itmp = $ftmp . '.isPcr';

  # make primer seq fasta file
  open(FPT, ">$ptmp");
  print FPT join("\t", 'PRIMER', $prm->forward_primer->seq->seq, $prm->reverse_primer->seq->seq);
  close(FPT);
  
  # run isPcr 
  system("./bin/isPcr", "$fa", "$ptmp", $outfmt, $itmp);

  # parse isPcr data
  open(FPCR, "<$itmp") or die "Cannot open isPcr temp file - $itmp!";
  # there should be only one line
  my $ln = <FPCR>;
  my @aln = split(/\t/, $ln);
  
  return ($aln[1], $aln[2]);
}

##### check primer compatibilty
# each primer trio will be grouped by unique common primer
sub primerSizeFloxCompatible {
  my ($fa, $prm0, $prm1) = @_;

  # get seq objects
  my $fain0 = $fa . '.0';
  my $seqio0 = Bio::SeqIO->new(-file => $fain0,
                               -format => 'fasta');
  my $seqobj0 = $seqio0->next_seq;
  my $fain1 = $fa . '.2';
  my $seqio1 = Bio::SeqIO->new(-file => $fain1,
                               -format => 'fasta');
  my $seqobj1 = $seqio1->next_seq;

  my $fain2 = $fa . '.1';
  my $seqio2 = Bio::SeqIO->new(-file => $fain2,
                               -format => 'fasta');
  my $seqobj2 = $seqio2->next_seq;  
  

  # parse 1st primer set to get unique common primers
  my %cmnSet = ();
  foreach my $pair0 (@$prm0) {
    my ($fp, $rp) = ($pair0->forward_primer, $pair0->reverse_primer);
    # check forward primer first
    if ($fp->oligo_type eq 'common') {
      my $seqstr = $fp->seq->seq;
      unless ( exists $cmnSet{$seqstr} ) {
        my @pps = ([], []);
        $cmnSet{$seqstr} = \@pps;
      }
      push @{$cmnSet{$seqstr}->[0]}, $pair0;
    }
    # then check reverse primer
    elsif ($rp->oligo_type eq 'common') {
      my $seqstr = $rp->seq->seq;
      unless ( exists $cmnSet{$seqstr} ) {
        my @pps = ([], []);
        $cmnSet{$seqstr} = \@pps;
      }
      push @{$cmnSet{$seqstr}->[0]}, $pair0;
    }
  }

  # parse 2nd primer set to group       
  foreach my $pair1 (@$prm1) {
    my ($fp, $rp) = ($pair1->forward_primer, $pair1->reverse_primer);
    my $fpseq = $fp->seq->seq;
    my $rpseq = $rp->seq->seq;

    # check forward primer first
    if (exists $cmnSet{$fpseq}) {
      push @{$cmnSet{$fpseq}->[1]}, $pair1;
    }
    elsif (exists $cmnSet{$rpseq}) {
      push @{$cmnSet{$rpseq}->[1]}, $pair1;
    }
  }

  # make combinatory primer sets and check multiplex scores
  my @prmResults = ();
  foreach my $sk (keys %cmnSet) {
    foreach my $ppr0 (@{$cmnSet{$sk}->[0]}) {
      foreach my $ppr1 (@{$cmnSet{$sk}->[1]}) {
        my @mscore = multiplxScore($ppr0, $seqobj0, $ppr1, $seqobj1);
        my $strn0 = shift(@mscore);
        my $strn1 = shift(@mscore);

        # prepare primer info in format:
        #   PrimPrimEnd2Score,PrimPrimEnd1Score,PrimPrimAnyScore
        #   CommPrimer,CommTm, WildPrimer/FloxPrimer,WildTm/FloxTm,
        #   WildProdSize,WildProdTm,FloxProdSize,FloxProdTm, 
        #   MutPrimer,MutTm,MutProdSize,MutProdTm
        #   WildProdSeq,FloxProdSeq,MutantProdSeq
        my @prmSet = convFloxSizePrimerSet($sk, $seqobj0, $ppr0, $strn0, $seqobj1, $ppr1, $strn1, $fain2, $seqobj2);
        push(@mscore, @prmSet);
        push @prmResults, \@mscore;
      }
    }
  }

  return @prmResults;
}

# primer input set format: seqObj_0, primerObj_0, strain_0, seqObj_1, primerObj_1, strain_1
# primer output set format: CommPrimer,CommTm, WildPrimer/FloxPrimer,WildTm/FloxTm, 
#                           WildProdSize,WildProdTm,FloxProdSize,FloxProdTm, 
#                           MutPrimer,MutTm,MutProdSize,MutProdTm
#                           WildProdSeq,FloxProdSeq,MutProdSeq
sub convFloxSizePrimerSet {
  my ($cmnSeq, $sqObj0, $pr0, $str0, $sqObj1, $pr1, $str1, $fas2, $sqObj2) = @_;
  my @result = ($cmnSeq);
  my ($fp0, $rp0) = ($pr0->forward_primer, $pr0->reverse_primer);
  my ($fp1, $rp1) = ($pr1->forward_primer, $pr1->reverse_primer);

  # run isPcr to get PCR fragment location in Flox mice sequence
  my ($flx_pst, $flx_pnd) = primerFloxIsPcr($fas2, $pr0);

  # check common primer
  # common primer is always in wild primer set: either fp0 or rp0
  if ($fp0->seq->seq eq $cmnSeq) {
    $result[1] = $fp0->melting_temp;
    # primer pair 0 gives wild fragment
    $result[2] = $rp0->seq->seq;
    $result[3] = $rp0->melting_temp;
    my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
    $result[4] = length($seqStr0);
    $result[5] = meltSim($seqStr0);
    # FloxProdSize,FloxProdTm (hxu 09/12/12)
    my $seqStr2 = $sqObj2->subseq($flx_pst, $flx_pnd);
    $result[6] = length($seqStr2);
    $result[7] = meltSim($seqStr2);
    # primer pair 1 gives mutant fragment
    $result[8] = $rp1->seq->seq;
    $result[9] = $rp1->melting_temp;
    my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
    $result[10] = length($seqStr1);
    $result[11] = meltSim($seqStr1);
    $result[12] = $seqStr0;
    $result[13] = $seqStr2; # FloxProdSeq (hxu 09/12/12)
    $result[14] = $seqStr1;
  }
  else {
    $result[1] = $rp0->melting_temp;
    # primer pair 0 gives wild fragment
    $result[2] = $fp0->seq->seq;
    $result[3] = $fp0->melting_temp;
    my $seqStr0 = $sqObj0->subseq($fp0->start, $rp0->end);
    $result[4] = length($seqStr0);
    $result[5] = meltSim($seqStr0);
    # FloxProdSize,FloxProdTm (hxu 09/12/12)
    my $seqStr2 = $sqObj2->subseq($flx_pst, $flx_pnd);
    $result[6] = length($seqStr2);
    $result[7] = meltSim($seqStr2);
    # primer pair 1 gives mutant fragment
    $result[8] = $fp1->seq->seq;
    $result[9] = $fp1->melting_temp;
    my $seqStr1 = $sqObj1->subseq($fp1->start, $rp1->end);
    $result[10] = length($seqStr1);
    $result[11] = meltSim($seqStr1);
    $result[12] = $seqStr0;
    $result[13] = $seqStr2; # FloxProdSeq (hxu 09/12/12)
    $result[14] = $seqStr1;
  }
  return @result;
}

 
# print out designed primers
sub prtFloxPrimers {
  my $hPrms = shift;
  my $mmNm = shift; # (hxu 07/12/2012: add $mmNm for mouse strain name)

  print FOT 'Mouse Genotyping Primers for Mouse Strain ~ ' . $mmNm, "\n\n";
  
  foreach my $k (sort keys %$hPrms) {
    my $prms = $hPrms->{$k};
    print FOT 'Primer set - ' . $k . ':', "\n";

    if ($k eq 'Size different primers') {
      print FOT join("\t", 'PrimPrimEnd2','PrimPrimEnd1','PrimPrimAny',
                           'CommPrimer','CommTm', 'WildPrimer/FloxPrimer','WildTm/FloxTm','WildProdSize','WildProdTm','FloxProdSize','FloxProdTm', 
                           'MutPrimer','MutTm','MutProdSize','MutProdTm', 'WildProdSeq','FloxProdSeq','MutantProdSeq'), "\n";
    }
    else {
      print FOT join("\t", 'PrimPrimEnd2','PrimPrimEnd1','PrimPrimAny',
                           'CommPrimer','CommTm', 'WildPrimer','WildTm','WildProdSize','WildProdTm', 
                           'MutPrimer','MutTm','MutProdSize','MutProdTm', 'WildProdSeq','MutantProdSeq'), "\n";
    }
    
    my @sorted = sort { $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] } @$prms;
    foreach my $pt (@sorted) {
      print FOT join("\t", @$pt), "\n";
    }
    print FOT "\n\n";
  }      
}


# (hxu 08/29/2012: add for flox mice)
# get coordinates of flox sites
sub getFloxsites {
  my ($inseq, $infrg, $inprm) = @_;
  my @myFlxs;

  ##### get the feature locations
  foreach my $myf (@$infrg) {
    my @mytgs = $myf->get_tag_values('name');
    # only one tag "name" for each fragment feature
    if ($mytgs[0] eq 'floxsite.1.before') {
      $myFlxs[0] = [$myf->start, $myf->end];
    }
    elsif ($mytgs[0] eq 'floxsite.1.flox') {
      $myFlxs[1] = [$myf->start, $myf->end];
    }
    elsif ($mytgs[0] eq 'floxsite.1.after') {
      $myFlxs[2] = [$myf->start, $myf->end];
    }
    elsif ($mytgs[0] eq 'floxsite.2.before') {
      $myFlxs[3] = [$myf->start, $myf->end];
    }
    elsif ($mytgs[0] eq 'floxsite.2.flox') {
      $myFlxs[4] = [$myf->start, $myf->end];
    }
    elsif ($mytgs[0] eq 'floxsite.2.after') {
      $myFlxs[5] = [$myf->start, $myf->end];
    }
  }

  ##### convert feature sets
  # (hxu 10/03/2012: fix the feature changes bugs)
  my $delength = $myFlxs[3][1] - $myFlxs[1][0] + 1;
  my @outfrg;
  foreach my $myf (@$infrg) {
    if ($myf->end < $myFlxs[1][0]) {
       # features before deletion start: floxsite.1.flox->start
       # do not change coordinates
       push @outfrg, $myf;
    }
    elsif ($myf->start > $myFlxs[3][1]) {
       # features after deletion end: floxsite.2.before->end
       # change coordinates
       my $tmpst = $myf->start - $delength;
       $myf->start($tmpst);
       my $tmpnd = $myf->end - $delength;
       $myf->end($tmpnd);
       push @outfrg, $myf;
    }
    # do nothing for the deleted features
    #  floxsite.1.flox, floxsite.1.after, change, floxsite.2.before
    #  $myFlxs[1][0] to $myFlxs[3][1]
  } 
  

  # construct sequence after flox deletion  from:
  #   left, floxsite.1.before, floxsite.2.flox, floxsite.2.after, right
  # the following sequences are deleted:
  #   floxsite.1.flox, floxsite.1.after, change, floxsite.2.before
  my $myseq = $inseq->subseq(1, $myFlxs[0]->[1]); # seqs: left, floxsite.1.before->end
  $myseq .= $inseq->subseq($myFlxs[4]->[0], $inseq->length); # seqs: floxsite.2.flox->start, floxsite.2.after, right
  my $outseq = Bio::Seq->new( -display_id => 'floxDelete',
                              -seq => $myseq );

  ##### update primer design option: TARGET
  #     and update primer design locations
  for (my $i = 0; $i < scalar(@$inprm); $i++) {
    if ($inprm->[$i] =~ m/TARGET=/) {
      my $spc = ($myFlxs[1][1] - $myFlxs[1][0] + 1)
              + ($myFlxs[4][1] - $myFlxs[4][0] + 1)
              + ($myFlxs[5][1] - $myFlxs[5][0] + 1);
      $inprm->[$i] = 'TARGET=' . $myFlxs[1][0] . ',' . $spc;
    }
    else {
      # other primer design parameter affected by flox deletion
      # that are INCLUDED, EXCLUDED, JUNCTION
      if ($inprm->[$i] =~ /JUNCTION=/) {
        my $tmpln = $inprm->[$i];
        $tmpln =~ s/JUNCTION=//;
        my @ajs = split(',', $tmpln);
        $inprm->[$i] = 'JUNCTION=';
        foreach my $jt (@ajs) {
          if ($jt < $myFlxs[1][0]) {
            $inprm->[$i] .= $jt . ',';
          }
          elsif ($jt > $myFlxs[3][1]) {
            $jt = $jt - $delength;
            $inprm->[$i] .= $jt . ',';
          }
        }
        # remove the last ','
        chop($inprm->[$i]) if (length($inprm->[$i]) > 9);
      }
      elsif ($inprm->[$i] =~ /(IN|EX)CLUDED=/) {
        my $tmpln = $inprm->[$i];
        $tmpln =~ s/(IN|EX)CLUDED=//;
        my @xcls = split(/\s+/, $tmpln);
        $inprm->[$i] =~ s/CLUDED=\.+/CLUDED=/;
        foreach my $xc (@xcls) {
          my ($strt,$lngt) = split(',', $xc);
          if ($strt < $myFlxs[1][0]) {
            if ( ($strt + $lngt) < $myFlxs[1][0] ) {
              $inprm->[$i] .= $xc . ' ';
            }
            elsif ( ($strt + $lngt) > $myFlxs[1][0] ) {
              if ( ($strt + $lngt) < $myFlxs[3][1] ) {
                $lngt = $myFlxs[1][0] - $strt;
                my $txc = $strt . ',' . $lngt;
                $inprm->[$i] .= $txc . ' ';
              }
              elsif ( ($strt + $lngt) >= $myFlxs[3][1] ) {
                $lngt = $lngt - $delength;
                my $txc = $strt . ',' . $lngt;
                $inprm->[$i] .= $txc . ' ';
              }
            }
          }
          elsif ( ($strt >= $myFlxs[1][0]) && ($strt < $myFlxs[3][1]) ) {
            if ( ($strt + $lngt) < $myFlxs[3][1] ) {
              # the entire region is in deleted sequence
              # do nothing 
            }
            elsif ( ($strt + $lngt) >= $myFlxs[3][1] ) {
              $strt = $myFlxs[3][1]; 
              $lngt = $strt + $lngt - $myFlxs[3][1];
              my $txc = $strt . ',' . $lngt;
              $inprm->[$i] .= $txc . ' ';
            }            
          }
          elsif ($strt >= $myFlxs[3][1]) {
            $inprm->[$i] .= $xc . ' ';
          }
        }
        # remove the last ' '
        chop($inprm->[$i]) if (length($inprm->[$i]) > 9);
      }
    }
  } 
    
  return ($outseq, \@outfrg, $inprm, \@myFlxs);
} 

################################################################################
##### Flox mice primer design <-- END
################################################################################
  


################################################################################
################################################################################
#######                                                                  #######
#######                       Input & Output files                       #######
#######                                                                  #######
################################################################################
################################################################################


=head1 Input File

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

[xuh2@ehsscilp01/jobs tmp]$ cat 5443252c4ee9452719104cbb288692c6.mfa.0
>Orai1gt|mutant ORDER=0; TARGET=887,4 8603,4; EXCLUDED=890,199 8405,199
AGTGTGGAAGAGGTAGGGATAGTGTTTGGTAGttcatttttgttggttttggtttttcca
gacaggggTGGTCTGTGTAGTCCTGGCTATCCTGGAACTTGCTTTATAGACCTGACTGGC
CTCGAACTCagagatcctcctgcctcttctcttgagtgctgggattaaaggcctgtgcca
ccacCAGTTGGTAGTtgaggtttttaaattttttaatgtgttttaaatttttcttggtcT
GGAGAGATAGTTAAGCAGGTAAAGATGCCTACTGCCAAAATTGACATtctgagtttgatc
cctaggACCcatgtggtggaaagaaagaaaagtgcCTCCCTgggcagaggcaggcagatc
TCTGAgctcaaggccagcctggtcCACAGAGCAAGTTCCAGGATGACCAGGGAAGGCTAC
ACAGAGAGACTCTGTCTCTcgaaaaagaaatttacaGCGTGCATAATATACCTAACTCTA
CCCGCTTTCACTATTATGATTGACAGTGGCAGTGTTCCCAATTGAAAGATATTTACTGTG
TACCAAGAACCcttatttgttaaaataattctgtgaATTAATGGGTCATACCTATTTTCT
ATATATGGTAAGGCTGGGAGACACTAACTTCCTAAGGACACACAGCTGATATAGATTCAc
GCTTGCTCTCCTCATCAATACTTTTTTGTGACTTCACAACCCAGAAAGGGTGGTGGGCTT
TGGCATTCCCAGAAATTGAGACTGAAGCCTGGCGTGGTGGCACgcgcctttaatcccagc
actcaggaggcagaggcaggcgGATTTCTGAGTTAgaggacagccagggctacacagaga
aaccctgtcttgaaaaacaaaaaacagaaaacaaaaaaaaagaaaaaagaaactgagaCT
GAGGCCACATGACCGATCTGATGGGAAGAGGACTAGGATTTGAGCGAGCCGttcactttc
ctcttcccctctcctAGGTAGCGATGGTGGAAGTCCAGCTGGACACAGACCATGACTACC
CACCAGGGTTGCTCATCGTCTTTAGTGCCTGCACCACAGTGCTAGTGGCCGTGCACCTGT
TTGCCCTCACAGAAAGAAGGTCTCAAGTTTTAGCCGGTAGCCCGGATGGCCTTTCCTGCA
GACCCCTACCACtttaccctttccctttgAAGgctttcccacaccaccctccacACTTGC
CCCAAACACTGCCAACTATGTAGGAGGAAGGGGTTGGGACTAACAGAAGAACCCGTTGTG
GGGAAGCTGTTGGGAGGGTCACTTTATGTTCTTGCCCAAGGTCAGTTGGGTGGCCTGCTT
CTGATGAGGTGGTCCCAAGGTCTGGGGTAGAAGGTGAGAGGGACAGGCCACCAAGgtcag
ccccccccccctaTCCCATAGGAGCCAGGTCCCTCTCCTGGACAGGAAGACTGAAGGGGA
GATGCCAGAGActcagtgaagcctggggTACCCTATTGGAGtccttcaaggaaacaaact
tggCCTCACCAggcctcagccttggctcCTCCTGGGAACTCTACTGCCCTTGGGATCCTA
CCGTTCGTATAGCATACATTATACGAAGTTATGTGATAGGCCTTTTAGCTACATCTGCCA
ATCCAtctcattttcacacacacacacaccactttccttctggtcAGTGGGCACATGTCC
AGCCTCAAGTTTATATCACCACCCCCAATGCCCAACACTTGTATGGCCTTGGGCGggtca
tccccccccccacccccagtatctgCAACCTCAAGCTAGCTTGGGTGCGTTGGTTGTGGA
TAAGTAGCTAGACTCCAGCAACCAGTAacctctgccctttctcctccaTGACAACCAGGT
CCCAGGTCCcgaaaaccaaagaagaagaACGCAGATCGCATCGATAACTTCGTATAGCAT
ACATTATACGAAGTTATCGCAGATCTGGACTCTAGAGGATCCCGTCGTTTTACAACGTCG
TGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGC
CAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCT
GAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCT
GGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGG
TTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGT
TCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCT
ACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTG
CAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAG
CGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGG
CAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTT
GCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGA
TTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCT
ACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTT
CGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAA
CGTCGAAAACCCGAAACTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGA
ACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGA
GGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGT
TAACCGTCACGAGCATCATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCA
GGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAA
CCATCCGCTGTGGTACACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAA
TATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACC
GGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGAT
CATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTG
GATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCAC
GGCCACCGATATTATTTGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGC
TGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGAT
CCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCA
GGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTC
GCTGattaaatatgatgaaaaCGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGA
TACGCCGAACGATCGCCAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCA
TCCAGCGCTGACGgaagcaaaacaccagcagCAGTTTTTCCAGTTCCGTTTATCCGGGCA
AACCATCGAAGTGACCAGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTG
GATGGTGGCGCTGGATGGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCC
ACAAGGTAAACAGTTGATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACT
CTGGCTCACAGTACGCGTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACAT
CAGCGCCTGGCAGCAGTGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTC
CCACGCCATCCCGCATCTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAA
GCGTTGGCAATTTAACCGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAA
ACAACTGCTGACGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGG
CGTAAGTGAAGCGACCCGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGG
CCATTACCAGGCCGAAGCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGT
GCTGATTACGACCGCTCACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAA
AACCTACCGGATTGATGGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAG
CGATACACCGCATCCGGCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCG
GGTAAACTGGCTCGGATTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTG
TTTTGACCGCTGGGATCTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGA
AAACGGTCTGCGCTGCGGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGA
CTTCCAGTTCAACATCAGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCA
TCTGCTGCACGCGGAAGAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGG
TGGCGACGACTCCTGGAGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTA
CCATTACCAGTTGGTCTGGTGTCAGGGGATCCCCCGGGCTGCAGCCAATATGGGATCGGC
CATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGG
CTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCAGC
GCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTGCA
GGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCT
CGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAGGA
TCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATGCG
GCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGCAT
CGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAAGA
GCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGACGG
CGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAATGG
CCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGACAT
AGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTCCT
CGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTTGA
CGAGTTCTTCTGAGCGGGACTCTGGGGTTCGAAATGACCGACCAAGCGACGCCCAACCTG
CCATCACGAGATTTCGATTCCACCGCCGCCTTCTATGAAAGGTTGGGCTTCGGAATCGTT
TTCCGGGACGCCGGCTGGATGATCCTCCAGCGCGGGGATCTCATGCTGGAGTTCTTCGCC
CACCCCCCGGATCTAAGCTCTAGATAAGTAATGATCATAATCAGCCATATCACATCTGTA
GAGGTTttacttgctttaaaaaaCCTCCCACACCTCCCCCTGAACctgaaacataaaatg
aatgcaattgttgttgttaacttgTTTATTGCAGCTTATAATGGttacaaataaagcaat
agcatCACAAATTTCACaaataaagcatttttttcaCTGCATTCTAGTTGTGGTTTGTCC
AAACTCATCAATGTATCTTATCATGTCTGGATCCGGGGGTACCGCGTCGAGAAGTTCCTA
TTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCGTCGACGCGAATTCGAATTCGT
AATCATGTCATAGctgtttcctgtgtgaaaTTGTTATCCGCTCACAATTCCACACAACAT
ACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATT
AATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTA
ATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTC
GCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAA
GGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAA
AGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCT
CCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGAC
AGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCC
GACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTC
TCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTG
TGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGA
GTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAG
CAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTA
CACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAG
AGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCggtggtttttttgtttg
caAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTAC
GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATC
AAAAAGGATCTTCACCTAGatccttttaaattaaaaatgaagttttaaatCAATCTAAAG
TATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTC
AGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTAC
GATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTC
ACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGG
TCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAG
TAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTC
ACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTAC
ATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAG
AAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTAC
TGTCATGCCATCCGTAAGAtgcttttctgtgactggTGAGTACTCAACCAAGTCATTCTG
AGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGC
GCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACT
CTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTG
ATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAaaaacaggaaggcaaaa
TGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTT
TCAAAAACTGAGACTGAGGCCACATGACCGATCTGATGGGAAGAGGACTAGGATTTGAGC
GAGCCGttcactttcctcttcccctctcctAGGTAGCGATGGTGGAAGTCCAGCTGGACA
CAGACCATGACTACCCACCAGGGTTGCTCATCGTCTTTAGTGCCTGCACCACAGTGCTAG
TGGCCGTGCACCTGTTTGCcCTCATGAtcagcacctgcatcctgCCCAACATCGAGGCTG
TGAGCAACGTCCACAACCTCAACTCGGTCAAAGAGTCACCCCACGAGCGCATGCATCGCC
ACATCGAgctggcctgggccttctCCACGGTCATCGGGACGCTGCTTTTCCTAGCAGAGG
TCGTGCTGCTCTGCTGGGTCAAGTTCTTACCTCTCAAGAGGCAAGCGGGACAGCCAAGCC
CCACCAAGCCTCCCgCTGAATCAGTCATCGTCGCCAACCACAGCGACAGCAGCGGCATCA
CCCCGGGTGAGGCGGCAGCCATTGCCTCCACCGCCATCAtggttccctgtggcctgGTTT
TTATCGTCTTTGCtGTTCACTTCTACCGcTCCCTGGTCAGCCATAAGACGGACCGGCAGT
TCCAGGAGCTCAATGAGCTGGCCGAGTTTGCCCGCTtgcaggaccagctggacCAC


=cut

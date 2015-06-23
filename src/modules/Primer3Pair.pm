#
# Perl module for handling Bio::Tools::Primer3Redux::PrimerPair problem
#   This module should wrap up "forward_primer" & "reverse_primer"
#   for using with Primer3Pair objects 
#
# Author: Hong Xu
#
# You may distribute this module under the same terms as perl itself

=head1 SYNOPSIS
use Primer3Pair;

# create and populate a gene project
$p3p = new Primer3Pair;
$p3p->forward_primer($fp);
$p3p->reverse_primer($rp);
......

=cut
 
package Primer3Pair;
 
use strict;
use vars qw($VERSION);

 
$VERSION = 1.0;
 
=head2 new

  Arg [1]    : none
  Example    : $p3p = Primer3Pair->new;
  Description: Constructor.  Creates a new Primer3Pair object
  Returntype : Primer3Pair
  Exceptions : none
  Caller     : general

=cut

sub new {
    my ($class) = @_;

    my $self = {
      'forward_primer' => ' ',
      'reverse_primer' => ' ',
    };
    bless $self,$class;

    return $self;
}


=head2 forward_primer

  Arg [1]    : (optional) reference to Bio::Tools::Primer3Redux::Primer
  Example    : my $fp = $p3p->forward_primer(); 
  Description: Getter/Setter for the forward_primer of the Bio::Tools::Primer3Redux::PrimerPair
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub forward_primer {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'forward_primer'} = $value;
    }
    return $self->{'forward_primer'};
}

=head2 reverse_primer

  Arg [1]    : (optional) reference to Bio::Tools::Primer3Redux::Primer
  Example    : my $rp = $p3p->reverse_primer(); 
  Description: Getter/Setter for the reverse_primer of the Bio::Tools::Primer3Redux::PrimerPair
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub reverse_primer {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'reverse_primer'} = $value;
    }
    return $self->{'reverse_primer'};
}

1;


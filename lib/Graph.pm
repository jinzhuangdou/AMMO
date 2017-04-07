# $Id$

# Copyright @ 2002 - 2010 The Institute for Genomic Research (TIGR).
# All rights reserved.
# 
# This software is provided "AS IS".  TIGR makes no warranties, express or
# implied, including no representation or warranty with respect to the
# performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  TIGR does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing the
# copyrights, patent rights or property rights of others.
# 
# This software program may not be sold, leased, transferred, exported or
# otherwise disclaimed to anyone, in whole or in part, without the prior
# written consent of TIGR.

package Graph;
{

=head1 NAME

TIGR::Foundation - TIGR Foundation object

=head1 SYNOPSIS

  use TIGR::Foundation;
  my $obj_instance = new TIGR::Foundation;

=head1 DESCRIPTION

This module defines a structure for Perl programs to utilize
logging, version reporting, and dependency checking in a simple way.

=cut

   BEGIN {
      require 5.006_00;                       # error if using Perl < v5.6.0  
   }

   #use strict;
   use Cwd;
   use Cwd 'chdir';
   use Cwd 'abs_path';
   use File::Basename;
   use File::Spec;
   use Getopt::Long;
   use IO::Handle;
   use POSIX qw(strftime);
   use Sys::Hostname;

   require Exporter;

   our @ISA;
   our @EXPORT;
   @ISA = ('Exporter');


   ## internal variables and identifiers
   our $REVISION = (qw$Revision$)[-1];
   our $VERSION = '1.41'; 
   our $VERSION_STRING = "$VERSION (Build $REVISION)";
   our @DEPEND = ();                          # there are no dependencies


   ## prototypes

   # Functional Class : general
   sub new($);
   sub addHit($$);
   sub getHit($$);


   

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut



   sub new($) {

      my $self = {};
      $self->{hit}    ={}; 
      $self->{hitmap} ={};
      $self->{hitNum} = 0;
      return $self;
   }


=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut


   sub addHit($$){

      my $self = shift ;
      my $hit  = shift ;
      my @aln  = split(/\s+|\t/, $hit);
      $self->{hitmap}{$aln[0]} = $hit ;
   };


=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut

   sub getMarkerNum($){

      my $self = shift ;
      my $inputID = shift;
      if (defined $inputID){
         return $self->{$inputID}->{markerNum} ;
      }
      else{
         return $self->{markerNum} ;
      }
      return undef; 
       
   };

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut


   sub getMarkerOfLg($){

      my $self = shift ;
      my $index= shift ;
      return $self->{$index} if (defined $self && defined $index);
      return undef; 

   };

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut


   sub getLgName($){

      my $self = shift ;
      my $index= shift ;
      return $self->{$index}->{name} if (defined $self && defined $index);
      return undef; 

   };

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut


   sub getMarkerPos($){

      my $self = shift ;
      my $index= shift ;
      return $self->{$index}->{pos} if (defined $self && defined $index);
      return undef; 

   };


=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub getLgOfMarker($) {

      my $self = shift ;
      my $index= shift ;
      my @markerID =() ;
      my $markerNum = $self->getMarkerNum($index) ;
      for(my $i=1;$i <= $markerNum; $i++){
          push @markerID, $self->{$index}->{$i}  ;
      }
      return @markerID if @markerID > 0 ;
      return undef; 

   }

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub getNextMarker($) {

      my $self = shift ;
      my $index= shift ;
      return $self->{$index}->{next} if (defined $self 
        && defined $index && $self->{$index}->{next} ne "NONE");
      return undef; 

   }

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub setCtgOfLg($) {
      my $self     = shift ;
      my $markerID = shift ;
      my $ctgID    = shift ;
      
      return 0 if (! defined $markerID || ! defined $ctgID);
      
      $self->{$markerID}->{ctgIndex} = $ctgID; 
  
      if (exists $self->{ctgIndex}->{$ctgID}){
         return 0;
      } 
      else{
        $self->{ctgIndex}->{$ctgID} ++ ;
        $self->{ctgNum} ++ ;
        $self->{ctgIndex}->{$self->{ctgNum}} =$ctgID ;
      }
      
   }

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub setDisOfLg($) {
      my $self        =shift ;
      my $markerID    =shift ;
      my $dis         =shift ;

      $lgID =$self->getLgName ($markerID);
    
      return 0 if (! defined $dis || ! defined $markerID || ! defined $lgID );
      $self->{dis} = $self->{dis} + $dis ;
      $self->{$lgID}->{dis} =  $self->{$lgID}->{dis} + $dis ;
      $self->{$markerID}->{dis} = $self->{$lgID}->{dis};
      return 1;
   }

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub getCtgOfLg($) {
      my $self     = shift ;
      my $markerID = shift ;
      return  $self->{$markerID}->{ctgIndex} if defined $self && defined $markerID ;
      return  undef;
     
   }

=item $obj_instance = new TIGR::Foundation;

This function creates a new instance of the TIGR::Foundation
object.  A reference pointing to the object is returned on success.  Otherwise,
this method returns undefined.

=cut 

   sub getDisOfLg($) {
      my $self     = shift ;
      my $markerID = shift ;
      return  $self->{$markerID}->{dis} if defined $self && defined $markerID ;
      return  undef;
     
   }




}
1;

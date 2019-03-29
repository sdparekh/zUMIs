package Approx;

=head1 NAME

String::Approx - match and substitute approximately (aka fuzzy matching)

=head1 SYNOPSIS

	use String::Approx qw(amatch asubstitute aregex);

=head1 DESCRIPTION

B<Approximate> is defined here as I<k-differences>.  One I<difference>
is an insertion, a deletion, or a substitution of one character.
The I<k> in the I<k-differences> is the maximum number of differences.

For example I<1-difference> means that a match is found if there is
one character too many (insertion) or one character missing (deletion)
or one character changed (substitution).  Those are I<exclusive or>s:
that is, I<not> one of each type of modification but I<exactly one>.

=head2 The default approximateness

The default approximateness is I<10 %> of the length of the
approximate pattern or I<at least 1>: I<0-differences> being the exact
matching which can be done very effectively using the usual Perl
function C<index()> or normal regular expression matching.

=head2 amatch

	use String::Approx qw(amatch);

	amatch("PATTERN");
	amatch("PATTERN", @LIST);
	amatch("PATTERN", [ @MODS ]);
	amatch("PATTERN", [ @MODS ], @LIST);

The PATTERN is B<a string>, not a regular expression.  The regular
expression metanotation (C<. ? * + {...,...} ( ) | [ ] ^ $ \w ...>)
will be understood as literal characters, that is, a C<*> means in
regex terms C<\*>, not I<"match 0 or more times">.

The LIST is the list of strings to match against the pattern.
If no LIST is given matches against C<$_>.

The MODS are the modifiers that tell how approximately to match.  See
below for more detailed explanation.  B<NOTE>: The syntax really is
C<[ @MODS ]>, the square brackets C<[ ]> must be in there and it is
B<not> a string, no quotes of any kind around the C<[ ]>.  It is
an anonymous array, see L<perlref>.  See below for MODS examples.

In scalar context C<amatch()> returns the number of successful
matches.  In list context C<amatch()> returns the strings that
had matches.

Example:

	use String::Approx qw(amatch);

	open(WORDS, '/usr/dict/words') or die;

	while (<WORDS>) {
	    print if amatch('perl');
	}

or the same ignoring case:

	use String::Approx qw(amatch);

	open(WORDS, '/usr/dict/words') or die;

	while (<WORDS>) {
	    print if amatch('perl', ['i']);
	}

=head2 asubstitute

	use String::Approx qw(asubstitute);

	asubstitute("PATTERN", "SUBSTITUTION");
	asubstitute("PATTERN", "SUBSTITUTION", @LIST);
	asubstitute("PATTERN", "SUBSTITUTION", [ @MODS ]);
	asubstitute("PATTERN", "SUBSTITUTION", [ @MODS ], @LIST);

The PATTERN is B<a string>, not a regular expression.  The regular
expression metanotation (C<. ? * + {...,...} ( ) | [ ] ^ $ \w ...>)
will be understood as literal characters, that is, a C<*> means in
regex terms C<\*>, not I<"match 0 or more times">.

Also the SUBSTITUTION is B<a string>, not a regular expression.  Well,
mostly.  I<Most of the> regular expression metanotation (C<.>, C<?>,
C<*>, C<+>, ...) will be not understood as literal characters, that
is, a C<*> means in regex terms C<\*>, not I<"match 0 or more times">.
The understood notations are

=over 8

=item	C<$`>

the part I<before> the approximate match

=item	C<$&>

the approximately matched part

=item	C<$'>

the part I<after> the approximate match

=back

The MODS are the modifiers that tell how approximately to match.  See
below for more detailed explanation.  B<NOTE>: Yes, the syntax is
really C<[ @MODS ]>, the square brackets C<[ ]> must be in there and
it is B<not> a string, no quotes of any kind around the C<[ ]>.  It
is an anonymous array, see L<perlref>.  See below for MODS examples.

The LIST is the list of strings to substitute against the pattern.
If no LIST is given substitutes against C<$_>.

In scalar context C<asubstitute()> returns the number of successful
substitutions.  In list context C<asubstitute()> returns the strings
that had substitutions.

Examples:

	use String::Approx qw(asubstitute);

	open(WORDS, '/usr/dict/words') or die;
	while (<WORDS>) {
	    print if asubstitute('perl', '($&)');
	}

or the same ignoring case:

	use String::Approx qw(asubstitute);

	open(WORDS, '/usr/dict/words') or die;
	while (<WORDS>) {
	    print if asubstitute('perl', '($&)', [ 'i' ]);
	}

=head2 aregex

	use String::Approx qw(aregex);

	aregex("PATTERN");
	aregex("PATTERN", [ @MODS ]);

B<Always returns a list>: the list of regular expressions that can be
used to approximately match the wanted pattern, given the possible
modifiers.  If the pattern is long or the modifiers call for large
approximateness, the list will have more than one element.  How to use
the element past the first one, is dependent on the intended
application of the regular expressions.

	my (@regex) = aregex('Obi-wan Kenobi');
	print $regex[0], "\n";
	print $regex[1], "\n";

If the pattern is to be passed onto an external utility such as
C<grep(1)> the B<std> modifier may be useful.  If that is applied, the
Perlisms C<(?:)> and <(?=)> are avoided.

	my (@regex) = aregex('Obi-wan Kenobi', [ 'std' ]);

=head2 Modifiers

The MODS argument in amatch(), asubstitute(), and aregex(), is an
anonymous array (see L<perlref>) of strings that control the matching
of PATTERN.  The first three modifiers, B<i>, B<g>, and B<?>, are the
usual regular expression match/substitute modifiers, the rest are
special for approximate matching/substitution.

=over 8

=item	i

Match/Substitute ignoring case, case-insensitively.

=item	g

Substitute I<globally>, that is, all the approximate matches, not just
the first one.

=item	?

Stingy match: instead of matching greedily, as much as possibly, match
as little as possible.

=item	std

Used usually only with the C<aregex()> interface.  Uses the C<()>
instead of the C<(?:)> and does not use the C<?=> in the result
pattern.  Useful if the pattern is going to be fed to an external tool
that does not understand these Perl extensions.  Produces shorter
regular expressions -- but they might match slower.

=item	I<k>

The maximum number of differences.
For example 2.

=item	II<k>

The maximum number of insertions.
For example 'I2'.

=item	DI<k>

The maximum number of deletions.
For example 'D2'.

=item	SI<k>

The maximum number of substitutions.
For example 'S2'.

=item	I<k>%

The maximum relative number of differences.
For example '10%'.

=item	II<k>%

The maximum relative number of insertions.
For example 'I5%'.

=item	DI<k>%

The maximum relative number of deletions.
For example 'D5%'.

=item	SI<k>%

The maximum relative number of substitutions.
For example 'S5%'.

=back

I<The regular expression modifiers> C<o m s x> I<are> B<not supported>
because their definitions for approximate matching are less than clear.

The relative number of differences is relative to the length of the
PATTERN, rounded up: if, for example, the PATTERN is C<'bouillabaise'>
and the MODS is C<['20%']> the I<k> becomes I<3>.

If you want to B<disable> a particular kind of difference you need
to explicitly set it to zero: for example C<'D0'> allows no deletions.

In case of conflicting definitions the later ones silently override,
for example:

	[2, 'I3', 'I1']

equals

	['I1', 'D2', 'S2']

=head1 EXAMPLES

The following examples assume the following template:

	use String::Approx qw(amatch asubstitute);

	open(WORDS, "/usr/dict/words") or die;
	while (<WORDS>) {
		# <---
	}

and the following examples just replace the above 'C<# E<lt>--->' line.

=head2 Matching from the C<$_>

=over 8

=item Match 'perl' with one difference

	print if amatch('perl');

The I<one difference> is automatically the result in this case because
first the rule of the I<10 %> of the length of the pattern ('C<perl>')
is used and then the I<at least 1> rule.

=item Match 'perl' with case ignored

	print if amatch('perl', [ 'i' ]);

The case is ignored in matching (C<i>).
B<NOTE>: this option halves the speed.

=item Match 'perl' with one insertion

	print if amatch('perl', [ '0', 'I1' ]);

The I<one insertion> is easiest achieved with first disabling any
approximateness (C<0>) and then enabling one insertion (C<I1>).

=item Match 'perl' with zero deletions

	print if amatch('perl', [ 'D0' ]);

The I<zero deletion> is easily achieved with simply disabling any
deletions (C<D0>), the other types of differences, the insertions and
substitutions, are still enabled.

=item Stingy matching 

	print if amatch('perl', [ '?' ]);

Match stingily: as little is matched as possible, as opposed to the
default greedy matching, where as much is matched as possible.

=item Substitute 'perl' approximately with HTML emboldening

	print if asubstitute('perl', '<B>$&</B>', [ 'g' ]);

All (C<g>) of the approximately matching parts of the input are
surrounded by the C<HTML> emboldening markup.

=item Stingy substitution

	print if asubstitute('perl', '<B>$&</B>', [ '?' ]);

Substitution is now stingy: as little is substituted as possible,
as opposed to the default greedy substitution, where as much is
substituted as possible.  When stingy the 'B<$&>' naturally
tends to be shorter than when greedy -- and the 'B<&`>' and
the 'B<$'>' respectively longer.

=back

=head2 Matching from a list

The above examples match against the default variable B<$_>.
The rest of the examples show how the match from a list.
The template is now:

	use String::Approx qw(amatch asubstitute);

	open(WORDS, "/usr/dict/words") or die;
	@words = <words>;
	# <---

and the examples still go where the 'C<# E<lt>--->' line is.

=over 8

=item Match 'perl' with one difference from a list

	@matched = amatch('perl', @words);

The C<@matched> contains the elements of the C<@words> that matched
approximately.

=item Substitute 'perl' approximately with HTML emphasizing from a list

	@substituted = asubstitute('perl', '<EM>$&</EM>', [ 'g' ], @words);

The C<@substituted> contains B<with all> (C<g>) B<the substitutions>
the elements of the C<@words> that matched approximately.

=back

=head1 DIAGNOSTICS

=over 8

=item amatch: $_ is undefined: what are you matching against?

=item asubstitute: $_ is undefined: what are you matching against?

These happen when you have nothing in C<$_> and try to C<amatch()> or
C<asubstitute()>.  Perhaps you are using the Perl option C<-e> but you
did forget the Perl option C<-n>?

=item amatch: too long pattern.

This happens when the pattern is too long for matching.

When matching long patterns, C<String::Approx> attempts to partition
the match.  In other words, it tries to do the matching incrementally
in smaller parts.

If this fails the above message is shown.  Please try using shorter
match patterns.

See below for L<LIMITATIONS/Pattern length> for more detailed
explanation why this happens.

=item asubstitute: too long pattern.

This happens when the pattern is too long for substituting.

The partitioning scheme explained above that is used for matching long
patterns cannot, sadly enough, be used substituting.

Please try using shorter substitution patterns.

See below for L<LIMITATIONS/Pattern length> for more detailed
explanation why this happens.

=back

=head1 TIPS

=over 8

=item transposes

To match transposed letters (as in "trasnposed") use substitutions = 2.

=item case ignorance

Avoid this, the speed is halved.

=item stingy matching

Also this tends to be somewhat slower.

=back

=head1 VERSION

Version 2.7.

=head1 LIMITATIONS

=head2 Fixed pattern

The PATTERNs of C<amatch()> and C<asubstitute()> are fixed strings,
they are not regular expressions.  The I<SUBSTITUTION> of
C<asubstitute()> is a bit more flexible than that but not by much.

=head2 Pattern length

The used approximate matching algorithm is B<very aggressive>.  In
mathematical terms it is I<O(exp(n) * k**2)>. This means that
when the pattern length and/or the approximateness grows the
matching or substitution take much longer time and memory.

For C<amatch()> this can be avoided by I<partitioning> the pattern,
matching it in shorter subpatterns.  This makes matching a bit slower
and a bit more fuzzier, more approximate.  For C<asubstitute()> this
partitioning cannot be done, the absolute maximum for the substitution
pattern length is B<19> but sometimes, for example it the approximateness
is increased, even shorter patterns are too much.  When this happens,
you must use shorter patterns.

=head2 Speed

We are still now at release 2.3 about 100 times slower than B<agrep>.

If you do not know what C<agrep> is: it is a program like the UNIX
grep for searching text from within files but it knows, among other
things, how to do approximate matching.  B<NOTE>: all these speeds
were measured in one particular system using one particular set of
tests: your mileage will vary.

=head2 Incompatibilities with C<String::Approx> I<v1.*>

If you have been using regular expression modifiers (B<i>, B<g>) you
lose.  Sorry about that.  The syntax simply is not compatible.  I had
to choose between having C<amatch()> match and C<asubstitute()>
substitute elsewhere than just in $_ I<and> the old messy way of
having an unlimited number of modifiers.  The first need won.

B<There is a backward compability mode>, though, if you do not want to
change your C<amatch()> and C<asubstitute()> calls.  You B<have> to
change your C<use> line, however:

	use String::Approx qw(amatch compat1);

That is, you must add the C<compat1> symbol if you want to be
compatible with the C<String::Approx> version 1 call syntax.

=head1 AUTHOR

Jarkko Hietaniemi C<E<lt>jhi@iki.fiE<gt>>

=head1 ACKNOWLEDGEMENTS

	Alberto Fontaneda C<E<lt>alberfon@ctv.esE<gt>>
	Dmitrij Frishman C<E<lt>frishman@mips.biochem.mpg.deE<gt>>
	Lars Gregersen C<E<lt>lars.gregersen@private.dkE<gt>>
	Helmut Jarausch C<E<lt>jarausch@IGPM.Rwth-Aachen.DEE<gt>>
	Håkan Kjellerstrand C<E<lt>hakank@netch.seE<gt>>
	Slaven Rezic C<E<lt>eserte@cs.tu-berlin.deE<gt>>
	Nathan Torkington C<E<lt>gnat@frii.comE<gt>>

Alberto Fontaneda and Dmitrij Frishman found a bug in long patterns,
suggested a test, and tested the patch.

Lars Gregersen saw String::Approx 2.2 and Håkan Kjellerstrand's MakeRegex
and that moment he experienced a genuine Eureka.  The result: up to thirty
times faster String::Approx.  (MakeRegex was used for completely other
purposes).

Helmut Jarausch noticed that 2.3 asubst() failed its test case in 5.004_50+.

Slaven Rezic found a bug in 2.5 and supplied a test case and a fix for it.

Nathan Torkington is to blame for the new API of release 2. :-)

=cut

#require 5;

use strict;
local $^W = 1;

use vars qw($PACKAGE $VERSION $compat1
	    @ISA @EXPORT_OK
	    %P @aL @dL @Pl %Pp);

$PACKAGE = 'String::Approx';
$VERSION = 2.7;

$compat1 = 0;

require Exporter;

@ISA = qw(Exporter);

@EXPORT_OK = qw(amatch asubstitute aregex);

# Catch the 'compat1' tag.

sub import {
    my $this = shift;
    my (@list, $sym);
    foreach $sym (@_) { $sym eq 'compat1' ? $compat1 = 1 : push(@list, $sym) }
    local $Exporter::ExportLevel = 1; 
    Exporter::import($this, @list);
}

sub _estimate {
    my ($l, $m) = @_;
    my $p = 5 ** ($m + 2);

    # Trust me, I know what I am doing.
    (3 * $p * $l ** 2 + (8 - $p) * $l - $p) / 8;
    # That was before the prefix optimizer and the stinginess.
    # Currently this is too aggressive (meaning too hasty
    # partitioning), the aggressiveness is overestimated,
    # both the optimizer and the stinginess should bring it down.
    # One of these lifetimes I will derive the new correct formula.
}

sub _common_suffix_find {
    my @list = @_;
    my ($s, $t) = (0, -1);
    my ($l, $c, $n, $min);
    
    $min = length($list[0]);
    foreach (@list[1..$#list]) {
	$l = length;
	$min = $l if $l < $min;
    }

    while ( $s < $min ) {
	$c = substr($list[0], $t, 1);
	# Wimp out on potential metacharacters.
	last if $c =~ /[\.\)\]\?]$/;
	$n = 1;
	foreach (@list[1..$#list]) {
	    last unless substr($_, $t, 1) eq $c;
	    $n++;
	}
	last if $n < @list;
	$s++;
	$t--;
    }

    $t == -1 ? '' : substr($list[0], $t+1);
}

# The _common_prefix_optimize() borrows heavily for the
# MakeRegex module by Håkan Kjellerstrand <hakank@netch.se>
# that from a list of words computes a fast regular expression
# that matches those words.  Changed (?:...) added (as noted
# by Lars Gregersen <lars.gregersen@private.dk>), added more
# [] instead of (), added lookahead (?=), added simple suffix
# optimization (most of this code is dead code for String::Approx,
# however, and is removed from here).

sub _common_prefix_optimize {
  my ($std, $p, @list)=@_;
  
  return "$p@list" if @list == 1;
  
  my ( %hash, $prefix, @all );

  foreach ( sort @list ) {
    $prefix = substr( $_, 0, 1 );
    push ( @{ $hash{ $prefix } }, length( ) > 1 ? substr( $_, 1 ) : '' );
  }
  
  my $question = 0;

  foreach ( sort keys %hash ) { # Recurse.
    my $comm = _common_prefix_optimize( $std, $_, @{ $hash{ $_ } } );
    $question = 1 if $comm eq "";
    push( @all, $comm );
  }

  my $paren = "";
  
  if ( @all == 1 ) {
    $paren = "@all";
  } else {
    my $maxlen = 0;  
    my $any    = 0;
    my $len;
    foreach ( @all ) {
	$len = length;
	$maxlen = $len > $maxlen ? $len : $maxlen;
	$any++ if /^\./;
    }

    my $mark = $question ? "?" : "";
    my $join;
    my $cando_lookahead = 0;
    my $suffix = @all ? _common_suffix_find( @all ) : '';
    # If we have pure-\w non-empty suffixes.
    if ( @all > 1 && $maxlen > 1 && length $suffix && not $suffix =~ /\W/ ) {
	my $sufflen = length $suffix;
	foreach ( @all ) { # Remove the suffix.
	    substr( $_, -$sufflen ) = '';
	}
	# Very heavy manual optimization took place here.
	# Two conditional branches out of three were removed.
	# I may have been wrong.  If that is the case,
	# this will produce illegal regexps.  That's life.
	# The cases that were removed apply in the general
	# case of doing prefix optimization, but not now.
	@all = grep length, @all;
	$join = "@all?$suffix";
    } else {
	my ( @l1, @lp );
	if ( $maxlen > 1 ) {
	    @l1 = grep { length() == 1 } @all;
	    @lp = grep { length() >  1 } @all;
	    push( @lp, @l1 == 1 ? @l1 : "[" . join( "", @l1 ) . "]" ) if @l1;
	} else {
	    @lp = @all;
	}
	$join = join( ($maxlen == 1) ? "" : "|", @lp );
	$cando_lookahead = 1;
    }
    if ( length $join == 1 ) {
	$paren = "$join$mark";
    } else {
      if ( $maxlen == 1 ) {
	  $paren = $any ? ".$mark" : "[$join]$mark";
      } else {
	  $paren = $std ? "($join)$mark" : "(?:$join)$mark";
	  if ( $any == 0 && $cando_lookahead && not $std ) {
	      my %first;
	      foreach ( @all ) {
		  $first{ substr( $_,0,1 ) } = undef;
	      }
	      my $class = join( '', sort keys %first );
	      if ( length $class ) {
		  $class = "[$class]" if length $class > 1;
		  $paren = "(?=$class)$paren";
	      }
	  }
      }
    }  
  }
  
  return "$p$paren";
}

sub _compile {
    my ($pattern, $I, $D, $S, $greed, $std) = @_;
    my ($j, $p, %p, %q, $l, $k, $mxm);
    my @p = ();

    $mxm = $I;	# maximum modifier
    $mxm = $D if ($D > $mxm);
    $mxm = $S if ($S > $mxm);

    $l = length($pattern);

    # The estimated length of the resulting pattern must be less than 32767.

    my $est = _estimate($l, $mxm);

    if ($est > 32767) { # The magic limit of Perl.
	my ( $a, $b, $i );
	my ( $mp, $np );

	$np = int( log( $l ) );

	# compute and cache the partitions per length

	unless (defined $Pl[$l][$mxm]) {
	    my ($sp, $fp, $gp);

	    $np = 2 if ($np < 2);
	    $sp = int($l / $np);
	    $fp = $l - $np * $sp;
	    $gp = $sp + $fp;
	    $mp = int($mxm / $np);
	    $mp = 1 if ($mp < 1);

	    $est = _estimate($gp, $mp);

	    while ($est > 32767) {
		# same rule here as above about the length of the pattern.
		$sp--;
		$np = int($l / $sp);
		$fp = $l - $np * $sp;
		$gp = $sp + $fp;
		$mp = int($mxm / $np);
		$mp = 1 if ($mp < 1);
		$est = _estimate($gp, $mp);
	    }

	    ($a, $b) = (0, $sp + $fp);
	    push(@{$Pl[$l][$mxm]}, [$a, $b]);
	    $a += $fp;
	    $b  = $sp;
	    for ($i = 1; $i < $np; $i++) {
		$a += $sp;
		push(@{$Pl[$l][$mxm]}, [$a, $b]);
	    }
	} else {
	    $mp = int($mxm / $np);
	}

	my $pi = $I ? int($mp / $I + 0.9) : 0;
	my $pd = $D ? int($mp / $D + 0.9) : 0;
	my $ps = $S ? int($mp / $S + 0.9) : 0;

	# compute and cache the pattern partitions

	unless (defined $Pp{$pattern}[$mxm]) {
	    foreach $i (@{$Pl[$l][$mxm]}) {
		push(@{$Pp{$pattern}[$mxm]},
		     [substr($pattern, $$i[0], $$i[1]), $pi, $pd, $ps]);
	    }
	}

	@p = @{$Pp{$pattern}[$mxm]};
	
    } else {
	push(@p, [$pattern, $I, $D, $S]);
    }

    my $i0 = 1;		# The start index for the insertions.

    my $pp;		# The current partition.

    foreach $pp (@p) {	# The partition loop.

	%p = ();

	my ($i, $d, $s) = @$pp[1..4];	# The per-partition I, D, S.

	$pp = $$pp[0];			# The partition string itself.

	$p{$pp} = length($pp);

	while ($i or $d or $s) {

	    %q = ();
	
	    # the insertions

	    if ($i) {
		$i--;
		while (($p, $l) = each %p) {
		    my $lp1 = $l + 1;

		    for ($j = $i0; $j < $l; $j++) {
			$k = $p;
			substr($k, $j) = '.' . substr($k, $j);
			$q{$k} = $lp1;
		    }
		}

		# After the first partition we want one insertion
		# before every partition - at index 0.  $i0 was
		# initialized before the partition loop as 1 and
		# thus the first partition does not get the one insertion
		# in front of it.

		$i0 = 0;
	    }

	    # the deletions

	    if ($d) {
		$d--;
		while (($p, $l) = each %p) {
		    if ($l) {
			my $lm1 = $l - 1;

			for ($j = 0; $j < $l; $j++) {
			    $k = $p;
			    substr($k, $j) = substr($k, $j + 1);
			    $q{$k} = $lm1;
			}
		    }
		}
	    }

	    # the substitutions

	    if ($s) {
		$s--;
		while (($p, $l) = each %p) {
		    for ($j = 0; $j < $l; $j++) {
			$k = $p;
			substr($k, $j, 1) = '.';
			$q{$k} = $l;
		    }
		}
	    }

	    while (($k, $l) = each %q) { $p{$k} = $l }
	}

	# If stingy, clean leading and trailing ".".

	unless ( $greed ) {
	    my (%c, $c);
	    foreach $c (keys %p) {
		$c =~ s/^\.+//;
		$c =~ s/\.+$//;
		$c{$c} = undef;
	    }
	    %p = %c;

	    my @lit = sort { length($b) <=> length($a) } keys %p;
	    my ( $li, $lj, %sh );

	    for ( $li = 0; $li < @lit; $li++ ) {
		for ( $lj = $li + 1; $lj < @lit; $lj++ ) {
		    if ( length $lit[ $lj ] < length $lit[ $li ]
			 and
			 $lj =~ /$li/ ) {
			$sh{ $lj } = undef;
		    }
		}
	    }

	    # Be pre-5.004-compatible.
	    foreach $c ( keys %sh ) {
		delete $p{ $c };
	    }
	}

	# The pattern.

	push(@{$P{$pattern}[$I][$D][$S][$greed][$std]},
             _common_prefix_optimize($std, "", keys %p));
    }
}

sub _mods {
    my ($mods, $aI, $aD, $aS, $rI, $rD, $rS) = @_;
    my $remods = '';
    my $mod;

    foreach $mod (@$mods) {
	while ($mod ne '') {
	    if ($mod =~ s/^([IDS]?)(\d+)(%?)//) {
		if ($1 ne '') {
		    if ($3 ne '') {
			if    ($1 eq 'I') { $$rI = 0.01 * $2 }
			elsif ($1 eq 'D') { $$rD = 0.01 * $2 }
			else              { $$rS = 0.01 * $2 }
		    } else {
			if    ($1 eq 'I') { $$aI = $2 }
			elsif ($1 eq 'D') { $$aD = $2 }
			else              { $$aS = $2 }
		    }
		} else {
		    if ($3 ne '') {
			$$rI = $$rD = $$rS = 0.01 * $2;
		    } else {
			$$aI = $$aD = $$aS = $2;
		    }
		}
	    } elsif ($compat1 and $mod =~ s/^([igmsxo])//) {
		$remods .= $1;
	    } elsif ($mod =~ s/^([ig])//) {
		$remods .= $1;
	    } elsif ($mod =~ s/^\?//) {
		# Just accept it.
	    } elsif ($mod =~ s/^std//) {
		# Just accept it.
	    } else {
		die $PACKAGE, ": unknown modifier '$mod'\n";
	    }
	}
    }

    $remods ne '' ? $remods : undef;
}

sub _mids {
    my ($len, $aI, $aD, $aS, $rI, $rD, $rS) = @_;

    my $r = int(0.1 * $len + 0.9);

    if    (    defined $rI) { $aI = int($rI * $len) }
    elsif (not defined $aI) { $aI = $r }

    if    (    defined $rD) { $aD = int($rD * $len) }
    elsif (not defined $aD) { $aD = $r }

    if    (    defined $rS) { $aS = int($rS * $len) }
    elsif (not defined $aS) { $aS = $r }

    ($aI, $aD, $aS);
}

sub _amatch {
    my ($pattern, $list) = @_;

    my ($aI, $aD, $aS, $rI, $rD, $rS);

    my $len = length($pattern);
    my $remods;
    my $greed = 1;
    my $std   = 0;

    if (ref $$list[0] or $compat1) {
	my $mods;

	if ($compat1) {
	    $mods = [ $$list[0] ];
	    @$list = ();
	} else {
	    $mods = shift(@$list);
	}

	$greed = 0 if grep { /\?/  } @$mods;
	$std   = 1 if grep { /std/ } @$mods;

	$remods = _mods($mods, \$aI, \$aD, \$aS, \$rI, \$rD, \$rS);

	($aI, $aD, $aS) = _mids($len, $aI, $aD, $aS, $rI, $rD, $rS);
    } else {
	$dL[$len] = int(0.1 * $len + 0.9) unless $dL[$len];
	$aI = $aD = $aS = $dL[$len];
    }

    _compile($pattern, $aI, $aD, $aS, $greed, $std)
	unless ref $P{$pattern}[$aI][$aD][$aS][$greed][$std];

    ( $len, $aD, $remods, @{$P{$pattern}[$aI][$aD][$aS][$greed][$std]} );
}

sub amatch {
    my ( $pattern, @list ) = @_;

    my ( $len, $aD, $remods, @mpat ) = _amatch( $pattern, \@list );

    die "amatch: \$_ is undefined: what are you matching against?\n"
	if (not defined $_ and @list == 0);

    my $mpat;

    @list = grep { length() > $len - $aD } @list if @list and $aD;

    if (@mpat == 1) {

	$mpat = $mpat[0];

	$mpat = '(?' . $remods . ')' . $mpat if defined $remods;

	if (@list) {

	    # match against the @list

	    my @m = eval { grep /$mpat/, @list };
	    die "amatch: too long pattern.\n" if ($@ =~ /regexp too big/);
	    return @m;
	}

	# match against the $_

	my $matched;

	eval { $matched = /$mpat/ };
	die "amatch: too long pattern.\n" if ($@ =~ /regexp too big/);
	return ($_) if $matched;

    } else {

	if (@list) {

	    # match against the @list

	    my @pos;
	    my @bad;
	    my ($i, $bad);

	    foreach $mpat (@mpat) {
		if (@pos) {
		    my $s;
		    foreach $s (@list) {
			pos($s) = shift(@pos);
		    }
		} else {
		    @pos = ();
		}
		for ($i = $bad = 0; $i < @list; $i++) {
		    unless (defined $bad[$i]) {
			if (eval { $list[$i] =~ /$mpat/g }) {
			    die "amatch: too long pattern.\n"
				if ($@ =~ /regexp too big/);
			    $pos[$i] = pos($list[$i]);
			} else {
			    $bad[$i] = $bad++;
			    return () if $bad == @list;
			}
		    }
		}
	    }
	    
	    my @got = ();

	    for ($i = 0; $i < @list; $i++) {
		push(@got, $list[$i]) unless defined $bad[$i];
	    }

	    return @got;
	}
	
	# match against the $_

	while ($mpat = shift(@mpat)) {
	    return () unless eval { /$mpat/g };
	    die "amatch: too long pattern.\n" if ($@ =~ /regexp too big/);
	    return ($_) if (@mpat == 0);
	}
    }

    return ();
}

sub _subst {
    my ($sub, $a, $b, $c) = @_;

    $sub =~ s/\$`/$a/g;
    $sub =~ s/\$&/$b/g;
    $sub =~ s/\$'/$c/g;

    $sub;
}

sub asubstitute {
    my ($pattern, $sub, @list) = @_;
    my ($aI, $aD, $aS, $rI, $rD, $rS);

    my $len = length($pattern);

    my $remods;
    my $greed = 1;
    my $std    = 0;

    if ($compat1 or ref $list[0]) {
	my $mods;

	if ($compat1) {
	    $mods = [ $list[0] ];
	    @list = ();
	} else {
	    $mods = shift(@list);
	}

	$greed = 0 if grep { /\?/  } @$mods;
	$std   = 1 if grep { /std/ } @$mods;

	$remods = _mods($mods, \$aI, \$aD, \$aS, \$rI, \$rD, \$rS);

	($aI, $aD, $aS) = _mids($len, $aI, $aD, $aS, $rI, $rD, $rS);
    } else {
	$dL[$len] = $len < 11 ? 1 : int(0.1 * $len) unless $dL[$len];
	$aI = $aD = $aS = $dL[$len];
    }

    die "asubstitute: \$_ is undefined: what are you matching against?\n"
	if (not defined $_ and @list == 0);

    _compile($pattern, $aI, $aD, $aS, $greed, $std)
	unless defined $P{$pattern}[$aI][$aD][$aS][$greed][$std];

    my @spat = @{$P{$pattern}[$aI][$aD][$aS][$greed][$std]};
    my $spat = $spat[0];
    
    $spat = '(?' . $remods . ')' . $spat if defined $remods;

    my ($s, $ss) = (0, 0);

    eval { '' =~ /$spat/ };
    die "asubstitute: too long pattern, maximum pattern length 19.\n"
	if ($@ =~ /regexp too big/);

    if (@list) {
	my ($sm, @m);

	foreach $sm (@list) {
	    # (?g) does not work.
	    if ( defined $remods and $remods =~ /g/ ) {
		push(@m, $sm)
		    if $sm =~ s/($spat)/_subst($sub,$`,$&,$')/eg;
	    } else {
		push(@m, $sm)
		    if $sm =~ s/($spat)/_subst($sub,$`,$&,$')/e;
	    }
	}

	return @m;
    }

    die "asubstitute: \$_ is undefined: what are you matching against?\n"
	unless defined $_;

    # (?g) does not work.
    if ( defined $remods and $remods =~ /g/ ) {
	return ($_)
	    if s/($spat)/_subst($sub,$`,$&,$')/eg;
    } else {
	return ($_)
	    if s/($spat)/_subst($sub,$`,$&,$')/e;
    }

    ();
}

sub aregex {
    my ( $pattern, @list ) = @_;
    
    my @mpat = _amatch( $pattern, \@list );

    splice( @mpat, 0, 3 ); # Junk first three.

    @mpat;
}

1;

# eof

#!/usr/bin/perl
#

sub PrntChty {
    my ( $chatty, $message ) = @_;
    use vars qw ( %Task );
    if ( $chatty <= $Task{chatter} ) {
        print"$message";
    }
    return;
} # PrntChty

sub GetInputParam {
    my ( $inspec ) = @_;
    use vars qw ( %Task );

    if ( $Task{status} ) { return; }

    my ( $ind ) = index ( $inspec, '=' );
    my ( $inpar ) = substr($inspec,0,$ind-1);
    my ( $inval ) = substr($inspec,$ind+1);

    if ( $ind <= 0 || !$inpar || ( !$inval && $inval ne "0" )) {
        print"$Task{'stem'}: Error: Parsing input parameter: $inspec\n";
        print"$Task{'stem'}: Error: Please specify parameters with format: <parametername>=<parametervalue> without blanks\n";
        print"$Task{'stem'}: Error: Type 'fhelp $Task{name}' for more information on parameters\n";
        $Task{status} = 1;
	    return;
    }
    
    return ($inval);

} #GetInputParam




1;

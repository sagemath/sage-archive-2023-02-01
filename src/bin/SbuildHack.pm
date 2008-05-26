# sbuildhack module originally from Debathena.
package SbuildHack;

use Sbuild qw(binNMU_version);

sub new_binNMU_version {
    my $v = shift;
    my $binNMUver = shift;
    die("Wrong binNMUver!") unless ($binNMUver == 171717);
    my %tags = (
        'sarge' => '~debian3.1',
        'etch' => '~debian4.0',
        'lenny' => '~debian4.1',
        'breezy' => '~ubuntu5.10',
        'dapper' => '~ubuntu6.06',
        'edgy' => '~ubuntu6.10',
        'feisty' => '~ubuntu7.04',
        'gutsy' => '~ubuntu7.10',
        'hardy' => '~ubuntu8.04'
    );
    return "$v$tags{$main::distribution}";
};

{
    no warnings 'redefine';
    *Sbuild::binNMU_version = \&new_binNMU_version;
}

1;

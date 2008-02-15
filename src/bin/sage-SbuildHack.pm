# sbuildhack module originally from Debathena.
package sage-SbuildHack;

use Sbuild qw(binNMU_version);
use Sbuild::Chroot qw(begin_session);
use Fcntl qw(:flock);

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

*old_begin_session = \&Sbuild::Chroot::begin_session;

sub new_begin_session {
    open(APTLOCK, ">/tmp/repository-lock");
    flock(APTLOCK, LOCK_SH);
    my $r = old_begin_session(@_);
    if (!open(PIPE, Sbuild::Chroot::get_apt_command("$conf::apt_get", "-q update", "root", 1) . " 2>&1 |")) {
	print PLOG "Can't open pipe to apt-get: $!\n";
	return 0;
    }
    while(<PIPE>) {
	print PLOG;
	print STDERR;
    }
    close(PIPE);
    print PLOG "apt-get update failed\n" if $?;
    flock(APTLOCK, LOCK_UN);
    close(APTLOCK);
    return $r;
}

{
    no warnings 'redefine';
    *Sbuild::binNMU_version = \&new_binNMU_version;
    *Sbuild::Chroot::begin_session = \&new_begin_session;
}

1;

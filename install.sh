#!/usr/bin/env sh

CPAN_MIRROR=https://mirrors.ustc.edu.cn/CPAN/
NO_TEST=--notest

hash cpanm 2>/dev/null || {
    echo >&2 "Install cpanm";
    curl -L https://cpanmin.us |
        $HOME/share/Perl/bin/perl - -v --mirror-only --mirror $CPAN_MIRROR App::cpanminus
}

cpanm --mirror-only --mirror $CPAN_MIRROR $NO_TEST Math::BigFloat IO::Tee Math::BigFloat Getopt::Long
cpanm --mirror-only --mirror $CPAN_MIRROR $NO_TEST Bio::SearchIO Bio::Seq Bio::SeqIO
cpanm --mirror-only --mirror $CPAN_MIRROR $NO_TEST AlignDB::IntSpan Sort::Array
cpanm --mirror-only --mirror $CPAN_MIRROR $NO_TEST File::Basename Path::Tiny Data::Dumper List::Util Getopt::Long

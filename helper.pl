#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Find 'find';
use File::Basename 'basename';
use File::Glob 'bsd_glob';

sub read_file {
  my $f = shift;
  open my $fh, "<", $f or die "FATAL: read_rawfile() cannot open file '$f': $!";
  binmode $fh;
  return do { local $/; <$fh> };
}

sub write_file {
  my ($f, $data) = @_;
  die "FATAL: write_file() no data" unless defined $data;
  open my $fh, ">", $f or die "FATAL: write_file() cannot open file '$f': $!";
  binmode $fh;
  print $fh $data or die "FATAL: write_file() cannot write to '$f': $!";
  close $fh or die "FATAL: write_file() cannot close '$f': $!";
  return;
}

sub check_source {
  my @all_files = (bsd_glob("makefile*"), bsd_glob("*.sh"), bsd_glob("*.pl"));
  find({ wanted=>sub { push @all_files, $_ if -f $_ }, no_chdir=>1 }, qw/src demo/);

  my $fails = 0;
  for my $file (sort @all_files) {
    next unless $file =~ /\.(c|h|pl|py|sh)$/ || basename($file) =~ /^makefile/i;
    my $troubles = {};
    my $lineno = 1;
    my $content = read_file($file);
    push @{$troubles->{crlf_line_end}}, '?' if $content =~ /\r/;
    for my $l (split /\n/, $content) {
      push @{$troubles->{merge_conflict}},   $lineno if $l =~ /^(<<<<<<<|=======|>>>>>>>)([^<=>]|$)/;
      push @{$troubles->{trailing_space}},   $lineno if $l =~ / $/;
      push @{$troubles->{tab}},              $lineno if $l =~ /\t/ && basename($file) !~ /^makefile/i;
      push @{$troubles->{non_ascii_char}},   $lineno if $l =~ /[^[:ascii:]]/;
      push @{$troubles->{cpp_comment}},      $lineno if $file =~ /\.(c|h)$/ && ($l =~ /\s\/\// || $l =~ /\/\/\s/);
      # in ./src we prefer using XMEMCPY, XMALLOC, XFREE ...
      push @{$troubles->{unwanted_memcpy}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmemcpy\s*\(/;
      push @{$troubles->{unwanted_malloc}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmalloc\s*\(/;
      push @{$troubles->{unwanted_realloc}}, $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\brealloc\s*\(/;
      push @{$troubles->{unwanted_calloc}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bcalloc\s*\(/;
      push @{$troubles->{unwanted_free}},    $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bfree\s*\(/;
      push @{$troubles->{unwanted_memset}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmemset\s*\(/;
      push @{$troubles->{unwanted_memcpy}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmemcpy\s*\(/;
      push @{$troubles->{unwanted_memmove}}, $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmemmove\s*\(/;
      push @{$troubles->{unwanted_memcmp}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bmemcmp\s*\(/;
      push @{$troubles->{unwanted_strcmp}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bstrcmp\s*\(/;
      push @{$troubles->{unwanted_strcpy}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bstrcpy\s*\(/;
      push @{$troubles->{unwanted_strlen}},  $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bstrlen\s*\(/;
      push @{$troubles->{unwanted_strncpy}}, $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bstrncpy\s*\(/;
      push @{$troubles->{unwanted_clock}},   $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bclock\s*\(/;
      push @{$troubles->{unwanted_qsort}},   $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bqsort\s*\(/;
      push @{$troubles->{sizeof_no_brackets}}, $lineno if $file =~ /^src\/.*\.c$/ && $l =~ /\bsizeof\s*[^\(]/;
      if ($file =~ m|src/.*\.c$| &&
          $l =~ /^static(\s+[a-zA-Z0-9_]+)+\s++([^s][a-zA-Z0-9_]+)\s*\(/) {
        push @{$troubles->{staticfunc_name}}, "$2";
      }
      if ($file =~ m|src/.*\.[ch]$| && $l =~ /^\s*#\s*define\s+(_[A-Z_][a-zA-Z0-9_]*)\b/) {
        my $n = $1;
        push @{$troubles->{invalid_macro_name}}, "$lineno($n)";
      }
      $lineno++;
    }
    for my $k (sort keys %$troubles) {
      warn "[$k] $file line:" . join(",", @{$troubles->{$k}}) . "\n";
      $fails++;
    }
  }

  warn( $fails > 0 ? "check-source:    FAIL $fails\n" : "check-source:    PASS\n" );
  return $fails;
}

sub check_comments {
  my $fails = 0;
  my $first_comment = <<'MARKER';
/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
MARKER
  my @all_files;
  find({ wanted=> sub { push @all_files, $_ if $_ =~ /\.(c|h)$/ }, no_chdir=>1 }, 'demo', 'src');
  for my $f (@all_files) {
    my $txt = read_file($f);
    if ($txt !~ /^\Q$first_comment\E/s) {
      warn "[first_comment] $f\n";
      $fails++;
    }
  }
  warn( $fails > 0 ? "check-comments:  FAIL $fails\n" : "check-comments:  PASS\n" );
  return $fails;
}

sub prepare_variable {
  my ($varname, @list) = @_;
  my $output = "$varname=";
  my $len = length($output);
  foreach my $obj (sort @list) {
    $len = $len + length $obj;
    $obj =~ s/\*/\$/;
    if ($len > 100) {
      $output .= "\\\n";
      $len = length $obj;
    }
    $output .= $obj . ' ';
  }
  $output =~ s/ $//;
  return $output;
}

sub patch_file {
  my ($content, @variables) = @_;
  for my $v (@variables) {
    if ($v =~ /^([A-Z0-9_]+)\s*=.*$/si) {
      my $name = $1;
      $content =~ s/\n\Q$name\E\b.*?[^\\]\n/\n$v\n/s;
    }
    else {
      die "patch_file failed: " . substr($v, 0, 30) . "..";
    }
  }
  return $content;
}

sub version_from_tfm_h {
  my $h = read_file(shift);
  if ($h =~ /\n#define\s*TFM_VERSION_S\s*"v([0-9]+)\.([0-9]+)\.([0-9]+)(\S*)"/s) {
    return "VERSION_PC=$1.$2.$3", "VERSION_LT=1:1", "VERSION=$1.$2.$3$4", "PROJECT_NUMBER=$1.$2.$3$4";
  }
  else {
    die "#define TFM_VERSION_S not found in tfm.h";
  }
}

sub make_sources_cmake {
  my ($list, $pub_headers) = @_;
  my $output = "set(SOURCES\n";

  foreach my $obj (sort @$list) {
    $output .= $obj . "\n";
  }
  $output .= ")\n\n";

  if ($pub_headers eq "") {
    return $output;
  }

  $output .= "set(PUBLIC_HEADERS\n";

  foreach my $obj (sort @$pub_headers) {
    $output .= $obj . "\n";
  }

  $output .= ")\n\nset(PRIVATE_HEADERS src/headers/tfm_private.h)\n";
  $output .= "set_property(GLOBAL PROPERTY PUBLIC_HEADERS \$\{PUBLIC_HEADERS\}\)\n\n";

  return $output;
}

sub process_makefiles {
  my $write = shift;
  my $changed_count = 0;
  my @c = ();
  find({ no_chdir => 1, wanted => sub { push @c, $_ if -f $_ && $_ =~ /\.c$/ && $_ !~ /.gen\.c$/ } }, 'src');
  my @h = ();
  find({ no_chdir => 1, wanted => sub { push @h, $_ if -f $_ && $_ =~ /\.h$/ && $_ !~ /tfm_private.h$/ } }, 'src');
  my @all = ();
  find({ no_chdir => 1, wanted => sub { push @all, $_ if -f $_ && $_ =~ /\.(c|h)$/  } }, 'src');

  my @o = sort (map { my $x = $_; $x =~ s/\.c$/.o/; $x } @c);
  my $var_o = prepare_variable("OBJECTS", @o);
  my $var_h = prepare_variable("HEADERS_PUB", (sort @h));
  (my $var_obj = $var_o) =~ s/\.o\b/.obj/sg;


  my @ver_version = version_from_tfm_h("src/headers/tfm.h");

  # update OBJECTS + HEADERS in makefile*
  for my $m (qw/ makefile makefile.shared sources.cmake /) {
    my $old = read_file($m);
    my $new = $m eq 'sources.cmake' ? make_sources_cmake(\@c, \@h)
            : patch_file($old, $var_o, $var_h, @ver_version);

    if ($old ne $new) {
      write_file($m, $new) if $write;
      warn "changed: $m\n";
      $changed_count++;
    }
  }

  if ($write) {
    return 0; # no failures
  }
  else {
    warn( $changed_count > 0 ? "check-makefiles: FAIL $changed_count\n" : "check-makefiles: PASS\n" );
    return $changed_count;
  }
}

sub die_usage {
  die <<"MARKER";
usage: $0 -s   OR   $0 --check-source
       $0 -c   OR   $0 --check-descriptors
       $0 -d   OR   $0 --check-defines
       $0 -o   OR   $0 --check-comments
       $0 -m   OR   $0 --check-makefiles
       $0 -a   OR   $0 --check-all
       $0 -u   OR   $0 --update-makefiles
       $0 --fixupind crypt.ind
MARKER
}

GetOptions( "s|check-source"        => \my $check_source,
            "c|check-descriptors"   => \my $check_descriptors,
            "d|check-defines"       => \my $check_defines,
            "o|check-comments"      => \my $check_comments,
            "m|check-makefiles"     => \my $check_makefiles,
            "a|check-all"           => \my $check_all,
            "u|update-makefiles"    => \my $update_makefiles,
            "f|fixupind=s"          => \my $fixupind,
            "h|help"                => \my $help
          ) or die_usage;

if ($fixupind) {
  my $txt = read_file($fixupind);
  $txt =~ s/^([^\n]*\n)/$1\n\\addcontentsline{toc}{chapter}{Index}\n/s;
  write_file($fixupind, $txt);
  exit 0;
}

my $failure;
$failure ||= check_source()       if $check_all || $check_source;
$failure ||= check_defines()      if $check_all || $check_defines;
$failure ||= check_comments()     if $check_all || $check_comments;
$failure ||= process_makefiles(0) if $check_all || $check_makefiles;
$failure ||= process_makefiles(1) if $update_makefiles;

die_usage unless defined $failure;
exit $failure ? 1 : 0;

# ref:         $Format:%D$
# git commit:  $Format:%H$
# commit time: $Format:%ai$

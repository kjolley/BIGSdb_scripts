#!/usr/bin/env perl
#Synchronise rMLST scratch database users with registered users
#Written by Keith Jolley, 2020-2022.
#Version 20220421
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR               => '/etc/bigsdb',
	LIB_DIR                  => '/usr/local/lib',
	DBASE_CONFIG_DIR         => '/etc/bigsdb/dbases',
	SCRATCH_DB_CONFIG        => 'pubmlst_rmlst_isolates_scratch',
	RMLST_ISOLATES_DB_CONFIG => 'pubmlst_rmlst_isolates',
	QUOTA                    => 200,
	MAX_AGE                  => 30 * 6
};
#######End Local configuration###############################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
my $scratch = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => SCRATCH_DB_CONFIG
	}
);
my $public = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => RMLST_ISOLATES_DB_CONFIG
	}
);
main();
undef $scratch;
undef $public;

sub main {
	delete_old_records();
	sync_admins();
	delete_users();
	add_users();
}

sub delete_old_records {
	my $age_days = MAX_AGE;
	eval {
		$scratch->{'db'}->do("DELETE FROM isolates WHERE datestamp < now()-interval '$age_days days'");
		$scratch->{'db'}->do('DELETE FROM retired_isolates');
	};
	if (@$) {
		$scratch->{'db'}->rollback;
		die qq($@\n);
	}
	$scratch->{'db'}->commit;
	return;
}

sub sync_admins {
	my $public_admins = $public->{'datastore'}
	  ->run_query( 'SELECT * FROM users WHERE status=?', 'admin', { fetch => 'all_arrayref', slice => {} } );
	my $scratch_usernames =
	  $scratch->{'datastore'}->run_query( 'SELECT user_name FROM users', undef, { fetch => 'col_arrayref' } );
	my %scratch_usernames = map { $_ => 1 } @$scratch_usernames;
	my $scratch_admin_usernames =
	  $scratch->{'datastore'}
	  ->run_query( 'SELECT user_name FROM users WHERE status=?', 'admin', { fetch => 'col_arrayref' } );
	my %scratch_admin_usernames = map { $_ => 1 } @$scratch_admin_usernames;
	foreach my $admin (@$public_admins) {
		next if $scratch_admin_usernames{ $admin->{'user_name'} };
		if ( $scratch_usernames{ $admin->{'user_name'} } ) {
			say qq(Setting $admin->{'user_name'} as admin.);
			eval {
				$scratch->{'db'}
				  ->do( 'UPDATE users SET status=? WHERE user_name=?', undef, 'admin', $admin->{'user_name'} );
			};
			if ($@) {
				$scratch->{'db'}->rollback;
				die qq($@\n);
			}
		} else {
			my $user_id = next_user_id();
			my @fields  = qw(user_name first_name surname affiliation email user_db);
			local $" = q(,);
			say qq(Adding $admin->{'user_name'} as admin.);
			eval {
				$scratch->{'db'}->do(
					"INSERT INTO users (id,@fields,status,date_entered,datestamp,curator) "
					  . 'VALUES (?,?,?,?,?,?,?,?,?,?,?)',
					undef, $user_id, @{$admin}{@fields}, 'admin', 'now', 'now', 0
				);
			};
			if ($@) {
				$scratch->{'db'}->rollback;
				die qq($@\n);
			}
		}
	}
	$scratch->{'db'}->commit;
}

sub delete_users {
	my $scratch_usernames =
	  $scratch->{'datastore'}->run_query( 'SELECT user_name FROM users', undef, { fetch => 'col_arrayref' } );
	my $public_usernames =
	  $public->{'datastore'}->run_query( 'SELECT user_name FROM users', undef, { fetch => 'col_arrayref' } );
	my %public = map { $_ => 1 } @$public_usernames;
	foreach my $user_name (@$scratch_usernames) {
		next if $public{$user_name};
		say qq(Deleting user $user_name.);
		eval {
			$scratch->{'db'}
			  ->do( 'DELETE FROM projects WHERE curator IN (SELECT id FROM users WHERE user_name=?)',
				undef, $user_name );
			$scratch->{'db'}->do( 'DELETE FROM users WHERE user_name=?', undef, $user_name );
		};
		if ($@) {
			$scratch->{'db'}->rollback;
			die qq($@\n);
		}
	}
	$scratch->{'db'}->commit;
}

sub add_users {
	my $scratch_usernames =
	  $scratch->{'datastore'}->run_query( 'SELECT user_name FROM users', undef, { fetch => 'col_arrayref' } );
	my $public_users =
	  $public->{'datastore'}
	  ->run_query( 'SELECT * FROM users ORDER BY user_name', undef, { fetch => 'all_arrayref', slice => {} } );
	my %scratch = map { $_ => 1 } @$scratch_usernames;
	my @fields = qw(user_name first_name surname affiliation email user_db);
	foreach my $user (@$public_users) {
		next if $scratch{ $user->{'user_name'} };
		my $user_id = next_user_id();
		local $" = q(,);
		say qq(Adding $user->{'user_name'}.);
		eval {
			$scratch->{'db'}->do(
				"INSERT INTO users (id,@fields,status,date_entered,datestamp,curator) "
				  . 'VALUES (?,?,?,?,?,?,?,?,?,?,?)',
				undef, $user_id, @{$user}{@fields}, 'submitter', 'now', 'now', 0
			);
			$scratch->{'db'}
			  ->do( 'INSERT INTO user_limits (user_id,attribute,value,curator,datestamp) VALUES (?,?,?,?,?)',
				undef, $user_id, 'private_isolates', QUOTA, 0, 'now' );
			my @permissions =
			  qw(modify_isolates modify_sequences tag_sequences designate_alleles delete_all only_private);
			foreach my $permission (@permissions) {
				$scratch->{'db'}->do( 'INSERT INTO permissions (user_id,permission,curator,datestamp) VALUES (?,?,?,?)',
					undef, $user_id, $permission, 0, 'now' );
			}
		};
		if ($@) {
			$scratch->{'db'}->rollback;
			die qq($@\n);
		}
	}
	$scratch->{'db'}->commit;
}

sub next_user_id {

	#this will find next id except when id 1 is missing
	my $next = $scratch->{'datastore'}->run_query(
		"SELECT l.id + 1 AS start FROM users AS l LEFT OUTER JOIN users AS r ON l.id+1=r.id "
		  . 'WHERE r.id is null AND l.id > 0 ORDER BY l.id LIMIT 1',
		undef,
		{ cache => 'next_user_id' }
	);
	$next = 1 if !$next;
	return $next;
}

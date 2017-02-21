#Migration of data between BIGSdb databases
#Written by Keith Jolley
#Copyright (c) 2011-2017, University of Oxford
package BIGSdb_Scripts::Migrate;
use strict;
use warnings;
use 5.010;
use Error qw(:try);
use List::MoreUtils qw(uniq);
use base qw(BIGSdb::Offline::Script);

sub run_script {
	my ($self) = @_;
	die "No connection to source database (check logs).\n" if !defined $self->{'db'};
	$self->_initiate_db( $self->{'options'}->{'b'} );
	if ( $self->{'options'}->{'c'} ) {
		my @client_dbs = split /,/, $self->{'options'}->{'c'};
		foreach (@client_dbs) {
			$self->_initiate_db($_);
			die "Client database $_ should be an isolate database.\n"
			  if $self->{'system2'}->{$_}->{'dbtype'} ne 'isolates';
		}
	}
	$self->{'user_map'} = $self->_map_users;
	return;
}

sub _initiate_db {
	my ( $self, $name ) = @_;
	my $xml_handler = BIGSdb::Parser->new;
	my $parser      = XML::Parser::PerlSAX->new( Handler => $xml_handler );
	my $full_path   = "$self->{'dbase_config_dir'}/$name/config.xml";
	eval { $parser->parse( Source => { SystemId => $full_path } ); };
	if ($@) {
		$self->{'logger'}->fatal("Invalid XML description: $@") if $self->{'instance'} ne '';
		die "Invalid database $name\n";
	}
	my %system_hash = %{ $xml_handler->get_system_hash };
	$self->{'system2'}->{$name} = $xml_handler->get_system_hash;
	$self->{'system2'}->{$name}->{'host'}     ||= 'localhost';
	$self->{'system2'}->{$name}->{'port'}     ||= 5432;
	$self->{'system2'}->{$name}->{'user'}     ||= 'apache';
	$self->{'system2'}->{$name}->{'password'} ||= 'remote';
	$self->_db_connect($name);
	return;
}

sub _db_connect {
	my ( $self, $name ) = @_;
	my %att = (
		'dbase_name' => $self->{'system2'}->{$name}->{'db'},
		'host'       => $self->{'system2'}->{$name}->{'host'},
		'port'       => $self->{'system2'}->{$name}->{'port'},
		'user'       => $self->{'system2'}->{$name}->{'user'},
		'password'   => $self->{'system2'}->{$name}->{'password'},
		'writable'   => 1
	);
	try {
		$self->{'db2'}->{$name} = $self->{'dataConnector'}->get_connection( \%att );
	}
	catch BIGSdb::DatabaseConnectionException with {
		$self->{'logger'}->error("Can not connect to database '$self->{'system2'}->{$name}->{'db'}'");
	};
	return;
}

sub _map_users {

	#Map users based on first and surnames.
	#Assumes that there is only one site-wide database.
	my ($self) = @_;
	my $db1_users = $self->{'datastore'}
	  ->run_query( 'SELECT * FROM users ORDER BY id', undef, { fetch => 'all_arrayref', slice => {} } );
	my $sql =
	  $self->{'db2'}->{ $self->{'options'}->{'b'} }
	  ->prepare('SELECT id, user_name FROM users WHERE first_name=? AND surname=?');
	my $sql_is_common_site_user = $self->{'db2'}->{ $self->{'options'}->{'b'} }
	  ->prepare('SELECT id FROM users WHERE user_name=? AND user_db IS NOT NULL');
	my $map;
	foreach my $db1_user (@$db1_users) {
		eval { $sql->execute( $db1_user->{'first_name'}, $db1_user->{'surname'} ) };
		$self->{'logger'}->error($@) if $@;
		my $db2_user = $sql->fetchrow_hashref;
		if ( $db2_user->{'user_name'} ) {
			$map->{ $db1_user->{'id'} } = $db2_user->{'id'};
		} else {
			eval { $sql_is_common_site_user->execute( $db1_user->{'user_name'} ) };
			$self->{'logger'}->error($@) if $@;
			my ($common_user) = $sql_is_common_site_user->fetchrow_array;
			if ($common_user) {
				$map->{ $db1_user->{'id'} } = $common_user;
			}
		}
	}
	return $map;
}

sub get_db_types {
	my ($self) = @_;
	return ( $self->{'system'}->{'dbtype'}, $self->{'system2'}->{ $self->{'options'}->{'b'} }->{'dbtype'} );
}

sub locus_exists_in_destination {
	my ( $self, $locus ) = @_;
	my $exists = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM loci WHERE id=?)",
		$locus,
		{ db => $self->{'db2'}->{ $self->{'options'}->{'b'} }, cache => 'Migrate::locus_exisgts_in_destination' } );
	return $exists;
}

sub locus_exists_in_source {
	my ( $self, $locus ) = @_;
	my $exists = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM loci WHERE id=?)",
		$locus, { cache => 'Migrate::locus_exists_in_source' } );
	return $exists;
}

sub _get_missing_curator_in_locus_tables {
	my ( $self, $locus ) = @_;
	my %bad_users;
	my @list;
	foreach my $table (
		qw (loci locus_extended_attributes locus_aliases locus_descriptions locus_links locus_refs locus_curators))
	{
		my $curators = $self->{'datastore'}->run_query(
			"SELECT DISTINCT "
			  . ( $table eq 'locus_curators' ? 'curator_id' : 'curator' )
			  . " FROM $table WHERE "
			  . ( $table eq 'loci' ? 'id' : 'locus' ) . "=?",
			$locus,
			{ fetch => 'col_arrayref' }
		);
		push @list, @$curators;
	}
	foreach ( uniq(@list) ) {
		my $user_info = $self->{'datastore'}->get_user_info($_);
		$bad_users{$_} = "$user_info->{'surname'}, $user_info->{'first_name'}" if !defined $self->{'user_map'}->{$_};
	}
	return \%bad_users;
}

sub get_missing_allele_seq_users_in_destination {
	my ( $self, $locus ) = @_;
	my $senders =
	  $self->{'datastore'}
	  ->run_query( "SELECT DISTINCT sender FROM sequences WHERE locus=?", $locus, { fetch => 'col_arrayref' } );
	my $curators =
	  $self->{'datastore'}
	  ->run_query( "SELECT DISTINCT curator FROM sequences WHERE locus=?", $locus, { fetch => 'col_arrayref' } );
	my %bad_users;
	foreach ( uniq( @$senders, @$curators ) ) {
		my $user_info = $self->{'datastore'}->get_user_info($_);
		$bad_users{$_} = "$user_info->{'surname'}, $user_info->{'first_name'}" if !defined $self->{'user_map'}->{$_};
	}
	return \%bad_users;
}

sub get_missing_designation_users_in_destination {
	my ( $self, $isolate_id ) = @_;
	my $senders = $self->{'datastore'}->run_query( "SELECT DISTINCT sender FROM allele_designations WHERE isolate_id=?",
		$isolate_id, { fetch => 'col_arrayref' } );
	my $curators =
	  $self->{'datastore'}->run_query( "SELECT DISTINCT curator FROM allele_designations WHERE isolate_id=?",
		$isolate_id, { fetch => 'col_arrayref' } );
	my %bad_users;
	foreach ( uniq( @$senders, @$curators ) ) {
		my $user_info = $self->{'datastore'}->get_user_info($_);
		$bad_users{$_} = "$user_info->{'surname'}, $user_info->{'first_name'}" if !defined $self->{'user_map'}->{$_};
	}
	return \%bad_users;
}

sub is_locus_in_scheme {
	my ( $self, $locus ) = @_;
	my $in_scheme = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM scheme_members WHERE locus=?)",
		$locus, { cache => 'Migrate::is_locus_in_scheme' } );
	return $in_scheme;
}

sub is_locus_in_destination {
	my ( $self, $locus ) = @_;
	my $in_destination = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM loci WHERE id=?)",
		$locus, { db => $self->{'db2'}->{ $self->{'options'}->{'b'} }, cache => 'Migrate::is_locus_in_destination' } );
	return $in_destination;
}

sub copy_locus {
	my ( $self, $locus ) = @_;
	die "Locus $locus does not exist in database $self->{'options'}->{'a'}\n" if !$self->locus_exists_in_source($locus);
	die "Locus $locus already exists in database $self->{'options'}->{'b'}\n"
	  if $self->locus_exists_in_destination($locus);
	die "Missing users in database $self->{'options'}->{'b'}"
	  if keys %{ $self->get_missing_allele_seq_users_in_destination($locus) };
	die "Missing users in database $self->{'options'}->{'b'}"
	  if keys %{ $self->_get_missing_curator_in_locus_tables($locus) };
	die "Locus $locus is a member of a scheme.\n" if $self->is_locus_in_scheme($locus);
	local $" = ',';
	eval {
		foreach my $table (
			qw(loci locus_extended_attributes locus_aliases locus_descriptions locus_links locus_refs locus_curators))
		{
			my ( @fields, @placeholders, @copy_data );
			my $attr = $self->{'datastore'}->get_table_field_attributes($table);
			my $data =
			  $self->{'datastore'}
			  ->run_query( "SELECT * FROM $table WHERE " . ( $table eq 'loci' ? 'id' : 'locus' ) . "=?",
				$locus, { fetch => 'row_hashref' } );
			if ( ref $data eq 'HASH' ) {
				$data->{'curator'}    = $self->{'user_map'}->{ $data->{'curator'} }    if defined $data->{'curator'};
				$data->{'curator_id'} = $self->{'user_map'}->{ $data->{'curator_id'} } if defined $data->{'curator_id'};
				foreach (@$attr) {
					push @fields,       $_->{'name'};
					push @placeholders, '?';
					push @copy_data,    $data->{ $_->{'name'} };
				}
				my $qry = "INSERT INTO $table (@fields) VALUES (@placeholders)";
				my $sql = $self->{'db2'}->{ $self->{'options'}->{'b'} }->prepare($qry);
				$sql->execute(@copy_data);
			}
		}
	};
	if ($@) {
		$self->{'logger'}->error($@);
		print "Can't insert locus data.\n";
		$self->{'db2'}->{ $self->{'options'}->{'b'} }->rollback;
	}
	$self->{'db2'}->{ $self->{'options'}->{'b'} }->commit;
	return;
}

sub copy_alleles {
	my ( $self, $locus ) = @_;
	my ( $sql_get, $sql_put, $fields );
	my @tables = qw (sequences sequence_extended_attributes accession sequence_refs allele_flags);
	foreach my $table (@tables) {
		my $attr = $self->{'datastore'}->get_table_field_attributes($table);
		my @placeholders;
		foreach (@$attr) {
			push @{ $fields->{$table} }, $_->{'name'};
			push @placeholders, '?';
		}
		local $" = ',';
		$sql_get->{$table} = $self->{'db'}->prepare("SELECT @{$fields->{$table}} FROM $table WHERE locus=?");
		my $qry = "INSERT INTO $table (@{$fields->{$table}}) VALUES (@placeholders)";
		$sql_put->{$table} = $self->{'db2'}->{ $self->{'options'}->{'b'} }->prepare($qry);
	}
	eval {
		foreach my $table (@tables) {
			$sql_get->{$table}->execute($locus);
			while ( my $data = $sql_get->{$table}->fetchrow_hashref ) {
				foreach (qw(curator sender)) {
					$data->{$_} = $self->{'user_map'}->{ $data->{$_} } if defined $data->{$_};
				}
				my @copy_data;
				push @copy_data, $data->{$_} foreach @{ $fields->{$table} };
				$sql_put->{$table}->execute(@copy_data);
			}
		}
	};
	if ($@) {
		$self->{'logger'}->error($@);
		print "Can't insert allele sequences.\n";
		$self->{'db2'}->{ $self->{'options'}->{'b'} }->rollback;
	}
	$self->{'db2'}->{ $self->{'options'}->{'b'} }->commit;
	return;
}

sub update_locus_fields_in_clients {
	my ( $self, $locus ) = @_;
	return if !$self->{'options'}->{'c'};
	my @client_dbs = split /,/, $self->{'options'}->{'c'};
	foreach my $client (@client_dbs) {
		my $new_name   = $self->{'system2'}->{ $self->{'options'}->{'b'} }->{'db'};
		my $old_name   = $self->{'system'}->{'db'};
		my $sql        = $self->{'db2'}->{$client}->prepare("UPDATE loci SET dbase_name=? WHERE id=? AND dbase_name=?");
		my $old_config = $self->{'options'}->{'a'};
		my $new_config = $self->{'options'}->{'b'};
		my $sql2 =
		  $self->{'db2'}->{$client}->prepare( "UPDATE loci SET url = replace(url,'$old_config',"
			  . "'$new_config'), description_url = replace(description_url,'$old_config','$new_config') WHERE id=?" );
		eval {
			$sql->execute( $new_name, $locus, $old_name );
			$sql2->execute($locus);
		};
		if ($@) {
			$self->{'logger'}->error($@);
			$self->{'db2'}->{$client}->rollback;
			print "Can't update clients.\n";
			return;
		}
		$self->{'db2'}->{$client}->commit;
	}
	return;
}

sub delete_locus_from_source {
	my ( $self, $locus ) = @_;
	my $sql  = $self->{'db'}->prepare("DELETE FROM sequences WHERE locus=?");
	my $sql2 = $self->{'db'}->prepare("DELETE FROM loci WHERE id=?");
	eval {
		$sql->execute($locus);
		$sql2->execute($locus);
	};
	if ($@) {
		$self->{'logger'}->error($@);
		$self->{'db'}->rollback;
		print "Can't delete locus $locus from source\n";
		return;
	}
	$self->{'db'}->commit;
	return;
}

sub clone_locus {
	my ( $self, $locus ) = @_;
	foreach my $table (qw (loci locus_aliases)) {
		my $attr = $self->{'datastore'}->get_table_field_attributes($table);
		my ( @fields, @placeholders );
		foreach (@$attr) {
			push @fields,       $_->{'name'};
			push @placeholders, '?';
		}
		local $" = ',';
		my $sql_get =
		  $self->{'db'}->prepare( "SELECT @fields FROM $table WHERE " . ( $table eq 'loci' ? 'id' : 'locus' ) . "=?" );
		my $qry     = "INSERT INTO $table (@fields) VALUES (@placeholders)";
		my $sql_put = $self->{'db2'}->{ $self->{'options'}->{'b'} }->prepare($qry);
		eval {
			$sql_get->execute($locus);
			my $data = $sql_get->fetchrow_hashref;
			if ( defined $data->{'curator'} ) {
				$data->{'curator'} = $self->{'user_map'}->{ $data->{'curator'} };
				my @copy_data;
				push @copy_data, $data->{$_} foreach @fields;
				$sql_put->execute(@copy_data);
			}
		};
		if ($@) {
			$self->{'logger'}->error($@);
			print "Can't clone locus $locus.\n";
			$self->{'db2'}->{ $self->{'options'}->{'b'} }->rollback;
			return;
		}
	}
	$self->{'db2'}->{ $self->{'options'}->{'b'} }->commit;
	return;
}

sub isolate_exists_in_destination {
	my ( $self, $isolate_id ) = @_;
	my $exists = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM isolates WHERE id=?)",
		$isolate_id,
		{ db => $self->{'db2'}->{ $self->{'options'}->{'b'} }, cache => 'Migrate::isolate_exists_in_destination' } );
	return $exists;
}

sub isolate_exists_in_source {
	my ( $self, $isolate_id ) = @_;
	my $exists = $self->{'datastore'}->run_query( "SELECT EXISTS(SELECT * FROM $self->{'system'}->{'view'} WHERE id=?)",
		$isolate_id, { cache => 'Migrate::isolate_exists_in_source' } );
	return $exists;
}

sub get_next_id {
	my ( $self, $table, $last_one ) = @_;
	my $next;
	if ( defined $last_one ) {
		$next = $last_one + 1;
		if ( !$self->{'sql'}->{'id_exists'}->{$table} ) {
			$self->{'sql'}->{'id_exists'}->{$table} =
			  $self->{'db2'}->{ $self->{'options'}->{'b'} }->prepare("SELECT EXISTS(SELECT id FROM $table WHERE id=?)");
		}
		eval { $self->{'sql'}->{'id_exists'}->{$table}->execute($next) };
		if ($@) {
			print $@;
			return;
		}
		my ($exists) = $self->{'sql'}->{'id_exists'}->{$table}->fetchrow_array;
		return $next if !$exists;
	}
	if ( !$self->{'sql'}->{'next_id'}->{$table} ) {

		#this will find next id except when id 1 is missing
		my $qry =
"SELECT l.id + 1 AS start FROM $table AS l left outer join $table AS r on l.id+1=r.id where r.id is null ORDER BY l.id LIMIT 1";
		$self->{'sql'}->{'next_id'}->{$table} = $self->{'db2'}->{ $self->{'options'}->{'b'} }->prepare($qry);
	}
	eval { $self->{'sql'}->{'next_id'}->{$table}->execute };
	if ($@) {
		print $@;
		return;
	}
	($next) = $self->{'sql'}->{'next_id'}->{$table}->fetchrow_array;
	$next = 1 if !$next;
	return $next;
}
1;

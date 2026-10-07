#!/usr/bin/env python3

"""
Create a dated snapshot of a PubMLST/BIGSdb scheme.

Written by Keith Jolley.
Copyright (c) University of Oxford 2026.
License: GPL3.

The resulting ZIP archive contains:

    alleles/
        <locus>.fasta
        ...
    profiles.tsv

The snapshot includes alleles and profiles entered on or before the
specified date.

As a consistency check, every allele referenced by a retained profile
is checked against the corresponding dated FASTA file. If an allele is
missing, a warning is printed and the allele is restored from the
current FASTA data.

Note that PubMLST/Pasteur BIGSdb will only allow retrieval data submitted
prior to 2025 if not authenticated. If you do not use an authentication key
then your results may not contain all the records that you expect.

    Set BIGSDB_API_KEY in the environment, or use --api-key.

"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import OrderedDict
from pathlib import Path
from typing import Any


DEFAULT_BASE_URL = "https://rest.pubmlst.org"
USER_AGENT = "bigsdb-scheme-snapshot/1.0"


class APIError(RuntimeError):
    """Raised when a BIGSdb API request fails."""


class BIGSdbClient:
    """Minimal serial BIGSdb REST client."""

    def __init__(
        self,
        base_url: str,
        api_key: str | None = None,
        timeout: int = 120,
        retries: int = 3,
    ):
        self.base_url = base_url.rstrip("/")
        self.api_key = api_key
        self.timeout = timeout
        self.retries = retries

    def api_url(self, path: str, **params: str) -> str:
        url = f"{self.base_url}/{path.lstrip('/')}"

        if params:
            url += "?" + urllib.parse.urlencode(params)

        return url

    def request(
        self,
        url: str,
        accept: str = "*/*",
    ) -> bytes:
        headers = {
            "Accept": accept,
            "User-Agent": USER_AGENT,
        }

        if self.api_key:
            headers["X-API-Key"] = self.api_key

        request = urllib.request.Request(url, headers=headers)

        for attempt in range(self.retries + 1):
            try:
                with urllib.request.urlopen(
                    request,
                    timeout=self.timeout,
                ) as response:
                    return response.read()

            except urllib.error.HTTPError as exc:
                try:
                    detail = exc.read().decode(
                        "utf-8",
                        errors="replace",
                    )
                except Exception:
                    detail = ""

                if exc.code < 500 or attempt == self.retries:
                    message = f"HTTP {exc.code} retrieving {url}"
                    if detail:
                        message += f": {detail}"
                    raise APIError(message) from exc

            except urllib.error.URLError as exc:
                if attempt == self.retries:
                    raise APIError(f"Error retrieving {url}: {exc}") from exc

            # Only reached for retryable errors.
            time.sleep(2**attempt)

        raise APIError(f"Failed to retrieve {url}")

    def get_text(
        self,
        url: str,
        accept: str = "text/plain",
    ) -> str:
        return self.request(url, accept).decode(
            "utf-8",
            errors="replace",
        )

    def get_json(self, url: str) -> Any:
        data = self.request(
            url,
            "application/json",
        )

        try:
            return json.loads(data)
        except json.JSONDecodeError as exc:
            raise APIError(f"Invalid JSON returned by {url}") from exc


def parse_date(value: str) -> dt.date:
    try:
        return dt.date.fromisoformat(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid date '{value}'; expected YYYY-MM-DD"
        ) from exc


def quote_component(value: str) -> str:
    return urllib.parse.quote(value, safe="")


def last_url_component(url: str) -> str:
    path = urllib.parse.urlparse(url).path.rstrip("/")
    return urllib.parse.unquote(path.split("/")[-1])


def safe_filename(name: str) -> str:
    return re.sub(
        r'[<>:"/\\|?*\x00-\x1f]',
        "_",
        name,
    )


def parse_fasta(
    text: str,
    locus: str,
) -> OrderedDict[str, str]:
    """
    Parse a BIGSdb allele FASTA into an ordered mapping:

        allele_id -> complete FASTA record

    BIGSdb FASTA identifiers are normally <locus>_<allele_id>,
    e.g. >pgm_231. The locus prefix is stripped so that dictionary
    keys correspond to the allele designations used in profiles.
    """

    records: OrderedDict[str, str] = OrderedDict()

    header: str | None = None
    lines: list[str] = []

    def save_record() -> None:
        if header is None:
            return

        fasta_id = header[1:].split()[0]

        prefix = f"{locus}_"
        if fasta_id.startswith(prefix):
            allele_id = fasta_id[len(prefix) :]
        else:
            allele_id = fasta_id

        record = "\n".join(lines)
        if not record.endswith("\n"):
            record += "\n"

        records[allele_id] = record

    for raw_line in text.splitlines():
        if raw_line.startswith(">"):
            save_record()
            header = raw_line
            lines = [raw_line]
        elif header is not None:
            lines.append(raw_line)

    save_record()

    return records


def write_fasta(
    filename: Path,
    records: OrderedDict[str, str],
) -> None:
    with filename.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as handle:
        for record in records.values():
            handle.write(record)


def parse_tsv(text: str) -> tuple[list[str], list[list[str]]]:
    """
    Parse the API's profiles_csv response while retaining the column
    order exactly as supplied by BIGSdb.
    """

    if not text.strip():
        return [], []

    reader = csv.reader(
        io.StringIO(text),
        delimiter="\t",
    )

    rows = list(reader)

    if not rows:
        return [], []

    return rows[0], rows[1:]


def write_tsv(
    filename: Path,
    header: list[str],
    rows: list[list[str]],
) -> None:
    with filename.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(header)
        writer.writerows(rows)


def get_scheme(
    client: BIGSdbClient,
    database: str,
    scheme_id: int,
) -> dict[str, Any]:

    url = client.api_url(f"db/{quote_component(database)}" f"/schemes/{scheme_id}")

    scheme = client.get_json(url)

    if not scheme.get("has_primary_key_field"):
        raise APIError(f"Scheme {scheme_id} does not have a primary key field.")

    return scheme


def get_field_name(
    client: BIGSdbClient,
    url: str,
) -> str:
    data = client.get_json(url)
    return str(data["field"])


def get_scheme_structure(
    client: BIGSdbClient,
    scheme: dict[str, Any],
) -> tuple[str, list[str]]:
    """
    Return the primary-key field name and ordered list of loci.
    """

    pk_url = scheme.get("primary_key_field")

    if not pk_url:
        raise APIError("Scheme does not provide a primary_key_field.")

    pk_field = get_field_name(client, pk_url)

    loci = [last_url_component(url) for url in scheme.get("loci", [])]

    if not loci:
        raise APIError("Scheme does not contain any loci.")

    return pk_field, loci


def get_locus_snapshot(
    client: BIGSdbClient,
    database: str,
    locus: str,
    snapshot_date: dt.date | None,
) -> tuple[
    OrderedDict[str, str],
    OrderedDict[str, str],
    int,
]:
    path = (
        f"db/{quote_component(database)}"
        f"/loci/{quote_component(locus)}"
        f"/alleles_fasta"
    )

    current_url = client.api_url(path)

    current_text = client.get_text(
        current_url,
        "text/plain",
    )

    current = parse_fasta(current_text, locus)

    # No date means a snapshot of all current records.
    if snapshot_date is None:
        return current.copy(), current, 0

    later_url = client.api_url(
        path,
        added_after=snapshot_date.isoformat(),
    )

    later_text = client.get_text(
        later_url,
        "text/plain",
    )

    later = parse_fasta(later_text, locus)

    snapshot = OrderedDict(
        (allele_id, record)
        for allele_id, record in current.items()
        if allele_id not in later
    )

    return snapshot, current, len(later)


def get_profile_snapshot(
    client: BIGSdbClient,
    database: str,
    scheme_id: int,
    snapshot_date: dt.date | None,
    pk_field: str,
) -> tuple[
    list[str],
    list[list[str]],
    int,
]:
    path = f"db/{quote_component(database)}" f"/schemes/{scheme_id}" f"/profiles_csv"

    current_url = client.api_url(path)

    current_text = client.get_text(
        current_url,
        "text/tab-separated-values",
    )

    header, current_rows = parse_tsv(current_text)

    if not header:
        raise APIError("Profile download returned no header.")

    # No date means retain the complete current profile set.
    if snapshot_date is None:
        return header, current_rows, 0

    later_url = client.api_url(
        path,
        added_after=snapshot_date.isoformat(),
    )

    later_text = client.get_text(
        later_url,
        "text/tab-separated-values",
    )

    later_header, later_rows = parse_tsv(later_text)

    try:
        pk_index = header.index(pk_field)
    except ValueError as exc:
        raise APIError(
            f"Primary key field '{pk_field}' is not present " "in profiles_csv."
        ) from exc

    if later_header and later_header != header:
        raise APIError(
            "Column layout of filtered profiles_csv differs "
            "from the unfiltered download."
        )

    later_ids = {row[pk_index] for row in later_rows if len(row) > pk_index}

    snapshot_rows = [
        row
        for row in current_rows
        if len(row) > pk_index and row[pk_index] not in later_ids
    ]

    return header, snapshot_rows, len(later_ids)


def split_allele_designations(value: str) -> list[str]:
    """
    Return allele identifiers represented by a profile cell.

    Normal BIGSdb profiles contain one allele identifier per locus
    column. This also tolerates semicolon-separated values.
    """

    value = value.strip()

    if not value:
        return []

    return [item.strip() for item in value.split(";") if item.strip()]


def ensure_profile_alleles(
    header: list[str],
    rows: list[list[str]],
    loci: list[str],
    snapshot_alleles: dict[
        str,
        OrderedDict[str, str],
    ],
    current_alleles: dict[
        str,
        OrderedDict[str, str],
    ],
) -> int:
    """
    Ensure every allele referenced by a retained profile is included
    in its corresponding FASTA file.

    Missing alleles are restored from the current allele FASTA.

    Returns the number of alleles restored.
    """

    column_index = {name: index for index, name in enumerate(header)}

    missing_locus_columns = [locus for locus in loci if locus not in column_index]

    if missing_locus_columns:
        raise APIError(
            "The following scheme loci are absent from profiles.tsv: "
            + ", ".join(missing_locus_columns)
        )

    restored = 0

    for locus in loci:
        index = column_index[locus]

        required: set[str] = set()

        for row in rows:
            if len(row) <= index:
                continue

            required.update(split_allele_designations(row[index]))

        dated = snapshot_alleles[locus]
        current = current_alleles[locus]

        for allele_id in required:
            if allele_id in dated:
                continue

            print(
                "WARNING: "
                f"{locus} allele {allele_id} is referenced "
                f"by a profile in the snapshot but was not "
                f"in the dated allele set; including it.",
                file=sys.stderr,
            )

            if allele_id not in current:
                # This is more serious: the profile refers to an
                # allele that isn't even in the current bulk FASTA.
                print(
                    "WARNING: "
                    f"{locus} allele {allele_id} could not be "
                    f"found in the current allele FASTA and "
                    f"cannot be added.",
                    file=sys.stderr,
                )
                continue

            dated[allele_id] = current[allele_id]
            restored += 1

    return restored


def make_zip(
    source_dir: Path,
    output_zip: Path,
) -> None:

    with zipfile.ZipFile(
        output_zip,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:

        for path in sorted(source_dir.rglob("*")):
            if not path.is_file():
                continue

            archive.write(
                path,
                path.relative_to(source_dir),
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a dated archive of a BIGSdb scheme "
            "containing allele FASTA files and profiles.tsv."
        )
    )

    parser.add_argument(
        "--database",
        required=True,
        help=("BIGSdb database configuration, e.g. " "pubmlst_neisseria_seqdef"),
    )

    parser.add_argument(
        "--scheme-id",
        required=True,
        type=int,
        help="Scheme ID",
    )

    parser.add_argument(
        "--date",
        type=parse_date,
        default=None,
        help=(
            "Snapshot date in YYYY-MM-DD format. "
            "If omitted, include all current records."
        ),
    )

    parser.add_argument(
        "--api-key",
        default=os.environ.get("BIGSDB_API_KEY"),
        help=(
            "BIGSdb personal API key. If omitted, " "BIGSDB_API_KEY is used when set."
        ),
    )

    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help=("REST API base URL " f"(default: {DEFAULT_BASE_URL})"),
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output ZIP filename",
    )

    args = parser.parse_args()

    today = dt.date.today()

    if args.date is not None and args.date > dt.date.today():
        parser.error(f"Snapshot date {args.date} is in the future.")

    if args.output is None:
        snapshot_date = args.date or dt.date.today()

        output = Path(
            f"{args.database}"
            f"_scheme_{args.scheme_id}"
            f"_{snapshot_date.isoformat()}.zip"
        )
    else:
        output = args.output

    output = output.resolve()

    client = BIGSdbClient(
        args.base_url,
        api_key=args.api_key,
    )

    if not args.api_key:
        print(
            "WARNING: No BIGSdb personal API key has been provided. "
            "The archive will contain only data accessible without "
            "authentication and may not include all expected records.",
            file=sys.stderr,
        )
        print(
            "Set BIGSDB_API_KEY or use --api-key to authenticate.",
            file=sys.stderr,
        )
        print(file=sys.stderr)

    print(
        f"Database:      {args.database}",
        file=sys.stderr,
    )
    print(
        f"Scheme:        {args.scheme_id}",
        file=sys.stderr,
    )
    snapshot_date = args.date or dt.date.today()

    print(
        f"Snapshot date: {snapshot_date}",
        file=sys.stderr,
    )
    print(
        f"API:           {args.base_url}",
        file=sys.stderr,
    )
    print(file=sys.stderr)

    print(
        "Retrieving scheme metadata...",
        file=sys.stderr,
    )

    scheme = get_scheme(
        client,
        args.database,
        args.scheme_id,
    )

    pk_field, loci = get_scheme_structure(
        client,
        scheme,
    )

    description = scheme.get(
        "description",
        f"scheme {args.scheme_id}",
    )

    print(
        f"Scheme: {description}",
        file=sys.stderr,
    )
    print(
        f"Primary key: {pk_field}",
        file=sys.stderr,
    )
    print(
        f"Loci: {len(loci)}",
        file=sys.stderr,
    )
    print(file=sys.stderr)

    with tempfile.TemporaryDirectory(prefix="bigsdb_snapshot_") as temp_name:

        root = Path(temp_name)
        allele_dir = root / "alleles"
        allele_dir.mkdir()

        snapshot_alleles: dict[
            str,
            OrderedDict[str, str],
        ] = {}

        current_alleles: dict[
            str,
            OrderedDict[str, str],
        ] = {}

        print(
            "Retrieving alleles...",
            file=sys.stderr,
        )

        for number, locus in enumerate(
            loci,
            start=1,
        ):
            dated, current, removed = get_locus_snapshot(
                client,
                args.database,
                locus,
                args.date,
            )

            snapshot_alleles[locus] = dated
            current_alleles[locus] = current

            print(
                f"[{number:>{len(str(len(loci)))}}/"
                f"{len(loci)}] "
                f"{locus}: "
                f"{len(dated)} retained, "
                f"{removed} later",
                file=sys.stderr,
            )

        print(
            "\nRetrieving profiles...",
            file=sys.stderr,
        )

        header, profile_rows, profiles_removed = get_profile_snapshot(
            client,
            args.database,
            args.scheme_id,
            args.date,
            pk_field,
        )

        print(
            f"Profiles retained: {len(profile_rows)}",
            file=sys.stderr,
        )
        print(
            f"Profiles later than snapshot: " f"{profiles_removed}",
            file=sys.stderr,
        )

        print(
            "\nChecking profile allele references...",
            file=sys.stderr,
        )

        restored = ensure_profile_alleles(
            header,
            profile_rows,
            loci,
            snapshot_alleles,
            current_alleles,
        )

        if restored:
            print(
                f"Added {restored} allele(s) required by " f"retained profiles.",
                file=sys.stderr,
            )
        else:
            print(
                "All profile allele references are present.",
                file=sys.stderr,
            )

        print(
            "\nWriting snapshot...",
            file=sys.stderr,
        )

        total_alleles = 0

        for locus in loci:
            filename = allele_dir / f"{safe_filename(locus)}.fasta"

            write_fasta(
                filename,
                snapshot_alleles[locus],
            )

            total_alleles += len(snapshot_alleles[locus])

        profiles_file = root / "profiles.tsv"

        write_tsv(
            profiles_file,
            header,
            profile_rows,
        )

        metadata = {
            "database": args.database,
            "scheme_id": args.scheme_id,
            "scheme_description": description,
            "snapshot_date": (
                args.date.isoformat()
                if args.date is not None
                else dt.date.today().isoformat()
            ),
            "generated": dt.datetime.now(dt.timezone.utc).isoformat(),
            "api_base_url": args.base_url,
            "primary_key": pk_field,
            "locus_count": len(loci),
            "allele_count": total_alleles,
            "profile_count": len(profile_rows),
            "alleles_restored": restored,
            "authenticated": bool(args.api_key),
        }

        metadata_file = root / "metadata.json"

        with metadata_file.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                metadata,
                handle,
                indent=2,
                ensure_ascii=False,
            )
            handle.write("\n")

        make_zip(root, output)

    print(
        "\nSnapshot complete.",
        file=sys.stderr,
    )
    print(
        f"Alleles:  {total_alleles}",
        file=sys.stderr,
    )
    print(
        f"Profiles: {len(profile_rows)}",
        file=sys.stderr,
    )
    print(
        f"Archive:  {output}",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(
            "\nInterrupted.",
            file=sys.stderr,
        )
        sys.exit(130)
    except (APIError, OSError) as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        sys.exit(1)

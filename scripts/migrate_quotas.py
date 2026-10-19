#!/usr/bin/env python3
"""Call the admin quota endpoint once for every row of an input TSV.

Input TSV (with header row) columns:
    user_email  pipeline_name  amount  quota_source  [comments]

`comments` is optional per row; --comments supplies a default for rows without one.

Example:
    ./migrate_quotas.py quotas.tsv \
        --base-url http://localhost:8080 \
        --token "$(gcloud auth print-access-token)" \
        --comments "TSPS-1255 quota migration" \
        --dry-run
"""
import argparse
import csv
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

REQUIRED_COLUMNS = ["user_email", "pipeline_name", "amount", "quota_source"]

def read_int(string_value):
    return int(string_value.replace(",",""))

def post_quota(base_url, token, row, comments):
    email = urllib.parse.quote(row["user_email"], safe="")
    pipeline = urllib.parse.quote(row["pipeline_name"], safe="")
    url = f"{base_url.rstrip('/')}/api/admin/v3/quotas/{pipeline}/{email}"
    body = json.dumps(
        {
            "quotaSource": row["quota_source"],
            "amount": read_int(row["amount"]),
            "comments": comments,
        }
    ).encode()
    req = urllib.request.Request(
        url,
        data=body,
        method="POST",
        headers={
            "accept": "*/*",
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.status, resp.read().decode()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv", help="input TSV file")
    parser.add_argument("--base-url", default="http://localhost:8080")
    parser.add_argument("--token", help="bearer token (default: $TSPS_TOKEN)")
    parser.add_argument("--comments", default="", help="comments for rows with no value in the comments column (required if any row lacks one)")
    parser.add_argument("--dry-run", action="store_true", help="print requests without sending them")
    parser.add_argument("--stop-on-error", action="store_true", help="abort at the first failed row")
    args = parser.parse_args()

    token = args.token or os.environ.get("TSPS_TOKEN")
    if not token and not args.dry_run:
        sys.exit("error: provide --token or set TSPS_TOKEN")

    with open(args.tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"error: input TSV missing columns: {', '.join(missing)}")
        rows = [{k: (v or "").strip() for k, v in r.items() if k} for r in reader]

    # The endpoint requires non-empty comments, so check every row before sending anything.
    no_comments = [i for i, r in enumerate(rows, start=2) if not (r.get("comments") or args.comments)]
    if no_comments:
        sys.exit(
            "error: comments are required; pass --comments or fill the comments column "
            f"(missing on line(s): {', '.join(map(str, no_comments))})"
        )

    failures = []
    for i, row in enumerate(rows, start=2):  # line 1 is the header
        label = f"line {i}: {row['user_email']} {row['pipeline_name']} {row['quota_source']} {read_int(row['amount'])}"
        try:
            if not all(row[c] for c in REQUIRED_COLUMNS):
                raise ValueError("empty required field")
            comments = row.get("comments") or args.comments
            if args.dry_run:
                print(f"DRY RUN  {label}")
                continue
            status, _ = post_quota(args.base_url, token, row, comments)
            print(f"OK ({status})  {label}")
        except urllib.error.HTTPError as e:
            detail = e.read().decode(errors="replace")
            print(f"FAILED ({e.code})  {label}: {detail}", file=sys.stderr)
            failures.append(label)
        except Exception as e:  # noqa: BLE001 - report and continue
            print(f"FAILED  {label}: {e}", file=sys.stderr)
            failures.append(label)
        if failures and args.stop_on_error:
            break

    print(f"\n{len(rows)} rows, {len(failures)} failed")
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()

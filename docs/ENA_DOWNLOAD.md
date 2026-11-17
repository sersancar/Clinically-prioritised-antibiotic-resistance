# ENA downloads with curl

The launch host must have `curl` on PATH. The downloader uses Python's standard
library plus the curl command; it no longer transfers FASTQs with urllib.
Run `curl --version` in the same shell used to launch Nextflow.

ENA HTTP access is attempted first, followed by HTTPS and FTP. Each transfer
must match the MD5 reported by the ENA run crosswalk before it is accepted.
There is no checksum bypass or TLS-verification bypass.

Interrupted transfers leave a `.fastq.gz.partial` file. Subsequent attempts in
the same task directory use curl's `--continue-at -` to resume from its byte
length. Three attempts per transport are allowed, with short backoff delays.
Connection establishment is limited to 30 seconds; transfers below 1 KiB/s
for 120 seconds are stopped and retried. There is no fixed total-download time
limit because these files are large.

If the server refuses byte-range resumption, the partial is removed and the
next attempt starts from zero. A completed transfer with a bad MD5 is also
removed and retried from zero. Previously verified destination files in that
task directory are reused.

Partial-file resumption applies within the same task directory. Nextflow may
create a different work directory for a process retry or after changing the
script; this helper does not search other task directories for partial files.
Keep `work/` and the Nextflow cache to reuse unchanged completed stages.

Validation: six tests use actual curl calls against a local HTTP server,
covering complete transfer, verified reuse, interrupted transfer resumption,
existing partial resumption, range refusal, HTTP 403 fallback and MD5 rejection.
A full production ENA FASTQ has not been downloaded during these tests.

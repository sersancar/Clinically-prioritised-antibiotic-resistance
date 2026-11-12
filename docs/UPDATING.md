# Apply this replacement to the existing GitHub checkout

The downloadable archive contains the new repository tree and an updater at
its top level. The updater compares the old repository files against the public
master commit used to prepare the replacement, then copies new files and removes
only the obsolete tracked paths. It does not change .git or push commits.

From the extracted archive, first check your existing checkout:

```bash
python apply_update.py /absolute/path/existing/repository
```

If the check passes, apply and review the changes:

```bash
python apply_update.py /absolute/path/existing/repository --apply
cd /absolute/path/existing/repository
git diff --stat
git add -A
git commit -m "Implement complete ENA-to-results analysis workflow"
git push
```

The checkout must be clean and previous files must match the expected contents.
The updater stops before changing anything if it finds local changes, unexpected
replacement files or a mismatched old version. Additional committed files that
are not part of this update, such as a licence, are preserved. Alternatively,
use the archive's clean source tree as a new working copy and review it manually.
Old scripts are absent from the new tree; this update does not rewrite Git history.

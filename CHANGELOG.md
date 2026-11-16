# Changelog

## 1.0.1 (2026-11-16)

Correction of the checksum file. No script, result file or number changed.

- In 1.0.0, Git converted the line endings of 53 result files (24 CSV, 23 JSON, 6 TXT) from CRLF to LF when they were
  committed, while `SHA256SUMS.txt` had been computed on the original files, so those 53 entries did not verify on a clone or
  on the release zip. 1.0.1 stores every file byte for byte (`.gitattributes`: `* -text`) and restores the original line
  endings, so the files are again identical to the files that were checksummed. Every checksum of 1.0.0 except the one for
  `README.md` is unchanged and now verifies.
- Two numba cache files (`__pycache__/*.nbc`, `*.nbi`) that had been committed by mistake are removed and ignored.
- `SHA256SUMS.txt` is computed on the committed contents. It now also lists `.gitattributes`, `.gitignore` and this file;
  the entry for `README.md` changed because the README describes this correction.

## 1.0.0 (2026-11-15)

First release, paper revision 15. Its checksum file does not verify for 53 files (see 1.0.1); those files differ from the
checksummed originals only in line endings. Kept unchanged.

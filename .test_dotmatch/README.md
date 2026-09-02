# DotMatch counting fixture

This directory exercises `count_method: dotmatch` on a 12-guide, fixed-length
library. The published Bassik test library is mixed-length (17–25 nt) and
cannot be used with DotMatch without first selecting one window.

Create the FASTQs before lint or a workflow run:

```bash
python3 workflow/scripts/make_dotmatch_test_reads.py
```

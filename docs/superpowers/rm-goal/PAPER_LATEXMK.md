# Paper LaTeX build transcript

Date: 2026-07-22 UTC

Worktree: /home/edgarcosta/projects/genus2isogenies/private/tex/_worktrees/rm-goal/

Branch: rm-goal-paper

Baseline commit: a92780134b7dfe25c91c60c42de09cf843fa6647

## Required command

Working directory:

~~~text
/home/edgarcosta/projects/genus2isogenies/private/tex/_worktrees/rm-goal/
~~~

Command, run without modification:

~~~sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
~~~

## Result

PASS. The first invocation performed a full build and exited with status 0.
It produced main.pdf with 21 pages and 536027 bytes.

The PDF SHA-256 was
ede2e538aca4e8568f3bf428756e72900feb51be0bd96a8f146ad905c1967795.
The new GonzalezGuardiaRotger2005 entry was present in main.bbl.

The build retained two unresolved citation warnings, Gro72 and BLR90.
Both citations occur in baseline main.tex, and neither key exists in the
baseline bibliography. They are unrelated to this insertion and did not make
the required command fail.

## Retained output from the full build

The terminal ended with the following actual output:

~~~text
Output written on main.pdf (21 pages, 536027 bytes).
Transcript written on main.log.
Latexmk: Getting log file 'main.log'
Latexmk: Examining 'main.fls'
Latexmk: Examining 'main.log'
Latexmk: Found input bbl file 'main.bbl'
Latexmk: Log file says output to 'main.pdf'
Latexmk: Using biber to make bibliography file(s).
Latexmk: Bibliography file(s) from .bcf file:
  genus2classes.bib
Sources for biber
  genus2classes.bib
  main.bcf
Latexmk: Summary of warnings from last run of *latex:
  Latex failed to resolve 2 citation(s)
Latexmk: All targets (main.pdf) are up-to-date
Latexmk: ====Undefined refs and citations with line #s in .tex file:
  Citation 'Gro72' on page 13 undefined on input line 682
  Citation 'BLR90' on page 13 undefined on input line 682
~~~

Exit status: 0.

## Complete output from the immediate verification rerun

The same required command was run again after the successful full build. Its
complete actual output was:

~~~text
Rc files read:
  /etc/LatexMk
Latexmk: This is Latexmk, John Collins, 15 June 2025. Version 4.87.
Latexmk: Nothing to do for 'main.tex'.
Latexmk: All targets (main.pdf) are up-to-date
~~~

Exit status: 0.

The generated untracked main.toc was removed after verification. It is not
part of the paper diff; the worktree has only the two intended tracked edits.

## Final orchestrator rerun

After Wave 3 synthesis, the orchestrator ran the same required command once
more.  Because `main.toc` had deliberately been removed, `latexmk` performed
a full convergent rebuild.  It again exited `0`, produced 21 pages and 536027
bytes, resolved the added `GonzalezGuardiaRotger2005` citation, and reported
only the same baseline `Gro72` and `BLR90` warnings.  The rebuilt PDF had
SHA-256
`d928bb5e517343e364b4a8cdb04d2347f404ea04b5dc04d9187717c8b10e974d`;
PDF build metadata makes that hash different from the earlier successful
build and it is not treated as a deterministic-output contract.  The newly
generated untracked `main.toc` was removed again.  Paper-worktree status then
contained only `main.tex` and `genus2classes.bib`.

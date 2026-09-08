# People Underestimate How Much Dissimilar Others Want to Connect

Data, code, and materials for:

Cooney, G.\*, Boothby, E. J.\*, Vorauer, J. D., & Schweitzer, M. E. (in press). People underestimate how much dissimilar others want to connect. *Journal of Personality and Social Psychology*.

\*Shared first authorship. Correspondence: Gus Cooney (gus.cooney@dartmouth.edu) or Erica Boothby (erica.j.boothby@dartmouth.edu).

## Summary

People systematically underestimate the extent to which dissimilar others are interested in socially connecting, a phenomenon we call *dissimilarity pessimism*. Across six experiments in the lab, online, and in the field, we document this pessimism across many kinds of dissimilarity (age, race, personality, job type, sociocultural background), explore two associated processes (greater uncertainty about others' interest and negatively skewed metaperceptual thoughts), and show that these beliefs are associated with people's choice of whom to approach, a previously underappreciated antecedent of homophily.

## Studies

| Study | Setting | Manipulation | Files |
|---|---|---|---|
| 1 | Lab | Imagined similar vs. dissimilar partner (general) | `study1_*` |
| 2 | Online | Same- vs. different-race partner (avatars) | `study2_*` |
| 3a | Online | Same vs. different personality type (inkblot task) | `study3a_*` |
| 3b | Online | Same vs. different personality type, with partner choice | `study3b_*` |
| 4 | Lab | Personality type and racial group; mechanism measures | `study4_*` |
| 5 | Field | Real colleagues ("Mystery Coffee"); five organizational dimensions | `study5_*` |

## Files

Each study has the same set of files:

- `studyX_AsPredicted.pdf` — preregistration
- `studyX_Qualtrics.qsf` — Qualtrics survey (importable); Study 2 has one file per partner-avatar condition
- `studyX_Qualtrics_Word.docx` — readable export of the survey
- `studyX_final.csv` — data (one row per participant; exclusions are applied in the scripts, not in the file)
- `studyX_final_public.R` — analysis script reproducing all reported statistics and figures for that study

Variable names in the CSVs follow the survey items and are self-descriptive; the `sim_`/`dissim_` and `self_`/`other_` prefixes index the similar/dissimilar partner condition and the self-rating/partner-rating, respectively. See the corresponding `_Qualtrics_Word.docx` for exact item wording.

Figure 5 (individual-level prevalence across studies) is produced in the "AGGREGATE" section at the end of `study4_final_public.R`. It reads the `studyX_long_aggregate.csv` files, each of which is written by the corresponding study script (search for `long_aggregate`).

## License

Code: MIT. Data and materials: CC BY 4.0. If you use these materials, please cite the paper above.

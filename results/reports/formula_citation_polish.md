# Formula rendering and citation attribution check

The previous PDF rendered the intended sample mean accent as an acute-like mark. The authoring source now uses `\overline{R}_{s}`, which produces a native Word upper-bar element, and the denominator is `\operatorname{sd}_{i}(R_{si})` to make its argument explicit. This is a notation and rendering correction; the implemented sample-wise mean and population standard deviation are unchanged.

Reference [8] is Doebley et al., Nature Communications 13, 7475 (2022), https://doi.org/10.1038/s41467-022-35076-w . Griffin is their software framework, so the original wording was not a bibliographic mismatch. The Introduction now attributes the finding directly to Doebley et al. for clarity, retaining the same reference and claim. The article describes nucleosome protection profiling, cancer phenotypes and fragment-length-specific GC correction.

The manuscript differs from the preceding version only in these two specified replacements. No classifier fitting, feature extraction or statistical recomputation was performed. The immutable v0.2.0 GitHub tag and Zenodo record are preserved; subsequent editorial corrections belong to the current manuscript and repository main branch.

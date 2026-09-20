# SA Score provenance

This directory vendors the RDKit Contrib `SA_Score` implementation so that
MinervaChem does not depend on the layout of a user's RDKit installation.

The files `sascorer.py` and `fpscores.pkl.gz` were copied from the local RDKit
installation at:

```text
rdkit/Contrib/SA_Score/
```

The implementation is based on:

> Peter Ertl and Ansgar Schuffenhauer, “Estimation of Synthetic Accessibility
> Score of Drug-like Molecules based on Molecular Complexity and Fragment
> Contributions,” Journal of Cheminformatics 1:8 (2009).

The vendored code retains its original copyright and BSD license notice in
`sascorer.py`. The score data in `fpscores.pkl.gz` is required by that code at
runtime and is distributed alongside it.

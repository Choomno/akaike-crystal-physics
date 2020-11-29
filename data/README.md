# Data

These three plain-text, space-delimited files are the inputs used in 2020. Values
are reproduced from the publications below for reference. Please cite the original
sources, not this repository.

## `swaminathan.txt`: urea, low temperature

Variables are stored in **rows** (`layout = 0`):

| row | quantity | unit |
|---|---|---|
| 1 | temperature *T* | K |
| 2 | lattice parameter *a* | Å |
| 3 | lattice parameter *c* | Å |

7 temperatures, 12–173 K. Urea is tetragonal, so *b* = *a*.

Source: S. Swaminathan, B. M. Craven et al., neutron-diffraction study of urea,
*Acta Crystallographica B* (1984). Structures also deposited in the Crystallography
Open Database (e.g. [COD 9017316](https://www.crystallography.net/cod/9017316.html)).

## `julian.txt`: urea, 188–328 K

Variables are stored in **columns** (`layout = 1`):

| column | quantity | unit |
|---|---|---|
| 1 | temperature *T* | K |
| 2 | lattice parameter *a* | Å |
| 3 | lattice parameter *c* | Å |

15 temperatures in 10 K steps.

Source: high-resolution variable-temperature X-ray powder diffraction by Hammond,
Pencheva, Roberts, Mougin and Wilkinson, *J. Appl. Cryst.* (2005),
[laboratory notes](https://journals.iucr.org/j/issues/2005/06/00/ko5017/index.html),
as tabulated in M. M. Julian, *Foundations of Crystallography with Computer
Applications* (hence the file name). There is a visible step in *a* between 248 K
and 258 K, and *a* at 298 K is slightly below *a* at 288 K.

## `anthracene.txt`: anthracene under hydrostatic pressure

Variables are stored in **columns** (`layout = 1`):

| column | quantity | unit |
|---|---|---|
| 1 | pressure *P* | GPa |
| 2 | lattice parameter *a* | Å |
| 3 | lattice parameter *b* | Å |
| 4 | lattice parameter *c* | Å |
| 5 | monoclinic angle β | degrees |

18 pressures, 1.05–22.7 GPa.

Source: synchrotron X-ray diffraction under pressure in a diamond anvil cell by
M. Oehzelt and co-workers (Graz University of Technology), who published this
work in several papers.

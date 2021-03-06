# Scripts for MATCH

## asts.py
Using the match fake file, you may be interested to know the 90% completeness
fraction.

- Load the ast file:
`ast = asts.AST(filename)`
- Compute the completeness
`ast.completeness(interpolate=True, combined_filters=True)`
- Compute the completeness fraction:
`ast.get_completeness_fraction(0.9)`

- Make a completeness plot with the 90% completeness marked:
`ast.completeness_plot(comp_fracs=[0.9])`

- Make the typical magin-magout plots:
`ast.magdiff_plot()`

You could also run it from the command line:
`python -m match.scripts.asts.py -h`


## calcsfh_parallel.py
Makes a script to run calcsfh in parallel. This will soon be depreciated or
vastly updated since it only uses taskset to run in paralell. If you aren't on
a unix machine, it's just a list of commands.

For now, try `python -m match.scripts.calcsfh_parallel.py -h` to see if anything breaks.


## diagnostics.py
This should be a very flexible entry to run all diagnostics on all files with
expected extensions in a directory. The extensions are currently:

- `.cmd` (cmd file that MATCH.calcsfh automatically makes)
- `.sfh` (calcsfh output)
- `.zc` (zcombine output (same format as .sfh))
- `.zc.dat` (zcmerge output (same format as .sfh with a different header))
- `.param` (match input parameter file)
- `.match` (match 2 column photometry)
- `.scrn` (calcsfh console output)


## fileio.py
utilities for opening, reading, writing files. The most useful is the `*_fmt()`
named functions. Those are the default files formates and any dictionary with
the correct keys can fill those values. It's a simple way to automatically write
parameter files. An example can be seen in `match.scripts.match_params.match_param`
TO DO: Make this a folder and separate functions. For example, all templates should be in their own file or together in a file.

## graphics.py
utilities for visualizing match outputs.

## likelihood.py
stellar_prob is a work in progress, lots of dumb warnings are thrown due to order of operations.


With main or run from command line: a call to MATCH `stats`.


## match_param.py
semi-automatic way to make `calcsfh` pararameter files.
Click on a cmd to set mag limits, or use completeness fraction, also can make
basic exclude gates.
It's set to read fits tables or match photometery. If it's a fits table, the
filters passed must exist exactly in the fits file.

## ssp.py
Work in progress. Visualizations and statistics on `calcsfh` in `ssp` mode

## utils.py
In flux. I think will be obsolete soon.

## cmd.py
CMD: was utils.MatchCMD (to read `.cmd` files, an automatic product from `calcsfh`)

## sfh.py
SFH: was utils.MatchSFH (to read `calcsfh` terminal output, `zcombine` and `zcmerge` output files)

## bash
Here are scripts to run MATCH

### write_hybridMC_script.sh
Used to run many hybridMC calls on all sfh files in a directory

### theworks.sh
Run all MATCH programs after a calcsfh (with mcdata flag) run

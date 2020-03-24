# SchistoIndividual
An R package providing wrapper functions for SchistoPkg.jl, an individual based model for schistosomiasis transmission.



# Installation

First seee platform specific instructions for installing Julia [here](https://julialang.org/downloads/platform/).

Then install SchistoIndividual using

```
devtools::install_github('mattg3004/SchistoIndividual')
```


Then you must setup the link between R and Julia using

```
SchistoIndividual()
```

If R cannot find your Julia install you might need to provide the path like this.

```
SchistoIndividual('/Applications/Julia-1.3.app/Contents/Resources/julia/bin/')
```

The first time you run this it will need to install a number of Julia packages and so will be quite slow.



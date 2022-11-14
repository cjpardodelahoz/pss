# PSS
R package to compute PSS index of phylogenetic diversity and the associated metrics (Pardo-De la Hoz et al., 2022).

## Installation
Install in R using devtools:

```R
install.packages("devtools") # Install devtools if you don't have it
library(devtools) # Load devtools
devtools::install_github("cjpardodelahoz/pss") # Install pss package
library(pss) # Load pss package
```

## Usage

### Get help within R

You can get a detailed help menu with a description of each function, and the respective arguments:

```R
?kl.mpd
?pss.r
?pss.rc
?ses.mpd.rc
```

### Calcualte KL-weighted mean pairwise distance (klMPD)

Let's use the data from Stanko et al. (2002), an interaction network between small mammals and their ecto-parasitic fleas from Slovakia. Load it like this:

```R
data("mf.samp") # The interaction matrix. Fleas are in the rows and mammals in the columns
data("mf.phyl.rows")  # Phylogeny of the fleas
data("mf.phyl.cols")  # Phylogeny of the mammals
```
We can get phylogenetic distance matrices from the phylogenetic trees with ape:

```R
install.packages("ape") # Install ape if you don't have it
library(ape) # load ape
mf.dis.rows <- cophenetic.phylo(mf.phyl.rows) # Phylogenetic distance matrix of fleas
mf.dis.cols <- cophenetic.phylo(mf.phyl.cols) # Phylogenetic distance matrix of the mammals
```

Now, let's calculate klMPD (see equation 8 in Pardo-De la Hoz et al. [2021]) for the fleas (rows) using the interaction matrix, and the distance matrix of the mammals (columns):

```R
kl.mpd(mf.samp, mf.dis.cols)
```

Notice that we didn't specify a value for the "q" argument. Therefore, the availability of the species was estimated from the interaction frequencies. Stanko et al. (2002) also quantified the availability of the mammal species directly in the field. Let's load those data:

```R
data("mf.q")
```

Now, we can calculate klMPD again, but with the user-specified vector of relative availabilities of the mammals:

```R
kl.mpd(mf.samp, mf.dis.cols, mf.q)
````

As you can tell, in this particular case, both approaches produce very similar results, which means that the interaction frequencies are a good proxy for the availability of the mammal species.

### Calculate PSS

Let's now calculate the PSS index for the fleas. This time, we can specify the phylogeny of the mammals directly. First, we can do it without specifying the vector with the empirical availabilities:

```R
pss.r(mf.samp, mf.phyl.cols)
```

And now, we calculate PSS using the empirical estimates of availability:

```R
pss.r(mf.samp, mf.phyl.cols, mf.q)
```
In both cases, the most important values are under the "pss" column of the resulting dataframe. values close to zero indicate a random phylogenetic structure, negative values indicate clustering, and positive values indicate overdispersion. See Pardo-De la Hoz et al. (2021) for more details on how to interpret this index.

If, as in this case, you have an interaction matrix and phylogenetic trees that include all taxa in the rows and the columns, you can get the PSS index for all taxa in a single line using pss.rc:

```R
pss.rc(mf.samp, mf.phyl.rows, mf.phyl.cols)
```

You can also specify a vector of relative availabilities in that function:

```R
pss.rc(mf.samp, mf.phyl.rows, mf.phyl.cols, q.cols = mf.q)
```

## Citation
Pardo-De la Hoz, C.J., Medeiros, I.D., Gibert, J.P., Chagnon, P.L., Magain, N., Miadlikowska, J., and Lutzoni, F. (2022). Phylogenetic structure of specialization: A new approach that integrates partner availability and phylogenetic diversity to quantify biotic specialization in ecological networks. Ecology and Evolution, 12, e8649.

## Contact
cjpardodelahoz@gmail.com\
carlos.pardo.de.la.hoz@duke.edu

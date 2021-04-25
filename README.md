# PSS
R package to compute PSS index of phylogenetic diversity and the associated metrics.

## Installation
Install in R using devtools:

```R
devtools::install_github("cjpardodelahoz/pss")
library(pss)
```

## Usage

### Calcualte KL-weighted mean pairwise distance (klMPD)

Let's use the data from Stanko et al. (2012), an interaction network between small mammals and their ecto-parasitic fleas from Slovakia. Load it like this:

```R
data("mf.samp") # The interaction matrix. Fleas are in the rows and mammals in the columns
data("mf.phyl.rows")  # Phylogeny of the fleas
data("mf.phyl.cols")  # Phylogeny of the mammals
```

Now let's calculate klMPD (see equation 8 in Pardo-De la Hoz et al. [2021]) for the fleas:

```R
kl.mpd(mf.samp, cophenetic(mf.phyl.cols))
#> Amalaraeus     Amphipsylla   Ceratophyllus Ctenocephalides  Ctenophthalmus 
      0.1422367      17.0621259      78.7565140      40.5100602      47.8661571 
    Dasypsyllus    Doratopsylla Hystrichopsylla     Leptopsylla     Megabothris 
      1.0258447      32.8802667      14.1511247      63.3680184      47.7557453 
    Nosopsyllus    Palaeopsylla Peromyscopsylla   Rhadinopsylla 
     14.2447009      24.6862258       0.7701843       2.5260813
```

Notice that we didn't specify a value for the "q" argument. Therefore, the availability of the species was estimated from the interaction frequencies. Stanko et al. (2012) also quantified the availability of the mammal species directly in the field. Let's load those data:

```R
data("mf.q")
```

Now, we can calculate klMPD again, but with the user-specified vector of relative availabilities of the mammals:

```R
kl.mpd(mf.samp, cophenetic(mf.phyl.cols), mf.q)
````

As you can tell, in this particular case, both approaches produce very similar results, which means that the interaction frequencies are a good proxy for the availability of the mammal species.

### Calculate PSS

Let's now calculate the PSS index for the fleas. First, we can do it without specifying the vector with the empirical availabilities:

```R
pss.r(mf.samp, mf.phyl.cols)
```

## Citation
Pardo-De la Hoz, CJ, Medeiros, ID, Gibert, JP, Chagnon, P, Magain, N, Miadlikowska, J, Lutzoni, F, (2021). A new approach to integrate phylogenetic structure and partner availability to study biotic specialization in ecological networks. IN REVIEW.

## Contact
cjpardodelahoz@gmail.com\
carlos.pardo.de.la.hoz@duke.edu

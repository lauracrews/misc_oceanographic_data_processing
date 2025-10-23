This repo contains miscellaneous useful functions for oceangraphic data processing.

`SKQ_underway_data_access.m` - Creates matlab structs of 1-minute bin averaged underway data from R/V Sikuliaq, accessing data from the [Coriolix API](https://coriolix.sikuliaq.alaska.edu/data/download/binned/). User makes API urls on Coriolix to specify desired parameters and time range, then copies and pastes these urls into the script. Optionally also renames variables and makes a struct compatible with [SWIFT codes](https://github.com/SASlabgroup/SWIFT-codes). 

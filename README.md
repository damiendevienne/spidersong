# SPIDER SONG

## Description

- Folder `prepareNHX` cointains R script (`code.R`) to prepare the NHX tree file. It also performs Brownian Motion simulation at nodes. It takes as input handmade chuncks of the original Nexus file provided by Wayne. 

-`SpiderFinal.nhx` is the tree file in NHX format generated by the above-mentionned R script. It is the input file for the `TreeSonif.py` script

- `TreeSonif.py` is a python script that transforms an input NHX tree file onto a json file compatible with the sonificatopn modules developed by Mendel. to execute the script simply type: 

```py
./TreeSonif.py -i input_file_name -o output_file_name
```

The ouput file name is optional. If not specified, the output is written to the standard output. 

To get help, type `TreeSonif.py -h`. To know the version of the script type `TreeSonif.py -v`.


## Checklist

- [x] Branch length associated to parent node
- [x] NeoY associated to parent node
- [x] Put in one place (function) all the features that are (or not) to be included in the json
- [ ] Add verbosity for people to know what happens when script is running
- [ ] Make the script general (1): choose what features to include (through option)
- [ ] Make the script general (2): choose what features to assign to parental nodes (through option) 
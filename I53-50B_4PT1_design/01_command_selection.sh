#!/bin/bash
rosetta_scripts.default.linuxgccrelease \
    -parser:protocol seletor_full.xml \
    -in:file:s preprocess/6p6f.pentamer.renumbered.pdb \
    -out:path:all selection

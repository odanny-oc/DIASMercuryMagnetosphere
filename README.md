# MESSENGER at Mercury: Tracking orbits and encounters with magnetospheric boundaries

This repository details the progress made on the analysis of the magnetospheric boundary crossings made by the MESSENGER spacecraft at Mercury from 2011 to 2015.

## Usage
To get started clone the repository using the link

```
git clone https://github.com/odanny-oc/DIASMercuryMagnetosphere.git
```

It is recommended to create a virtual environment inside the repository to download the requirements. To do this cd into the repository and run

```
cd DIASMercuryMagnetosphere
python -m venv .venv
```
To access the Python virtual environment run

```
source .venv/bin/activate
```

All the Python dependencies can be installed from the requirements.txt file in the code directory. To properly install the custom packages, run

```
cd code 
pip install -r requirements.txt
```

Once completed, the code, which all sits in the code directory, can be run from the terminal using

```
python /path/to/file
```

## Overview
This repository is a collection of code that produces the analysis of the Hollman et al. crossing list 2026 (https://zenodo.org/records/21392216). A detailed README can be found within each directory explaining how each file works. The analysis largely pertains to crossings and encounters of _magnetospheric boundaries_ of Mercury. There are two magnetospheric boundaries at Mercury, the _magnetopause_ (MP), and the _bow shock_ (BS). The bow shock is defined as the boundary where particles in the _solar wind_ are suddenly decelerated after first coming into contact with the magnetic field of Mercury. The magnetopause, which sits beyond the bow shock, is where the _magnetic pressure_ of the Sun's magnetic field equalises with Mercury's. This gives rise to three regions, the _solar wind_, the _magnetosheath_, and the _magnetosphere_.

A crossing, which is taken from the Hollman et al. dataset, is the moment when the spacecraft passes from one region to another. These then come with four different labels, MP_IN, MP_OUT, BS_IN, and BS_OUT.

# Warping
This repository contains the code morphing for warping a source mesh by deforming a closed template mesh. 

### Installation
You can grab the sources by cloning this repository or downloading a .zip archive of the sources. In order to build the project, navigate to the created directory and in a command prompt type:

mkdir build

cd build

cmake ..

make

(sudo) make install

### Usage
In a terminal, run:

usage: warping [-h] [-p] [- nit n] [-load l] [-lame lambda mu] [-yp young poisson] [-t template_file[.mesh]] -s source_file[.mesh]

The square braces indicate optional arguments. Some commands have flags, some others do not.

The options and flags are:

-h                         show the default parameters and exit

-p                         print output at each iteration

-nit  n                    number of iterations desired

-load l                    pressure magnitude 

-lame lambda mu            Lame coefficients 

-yp   young poisson        Young and Poisson coefficients
  

-t template_file.mesh      name of the template mesh file ( optional )

-s source_file.mesh        name of the source mesh file ( mandatory )

When the -t template_file[.mesh] is not provided, the code takes in input a default spherical template "sphere.mesh" (provided in the folder demo). 

### Authors and contributors

warping has been initiated by Maya de Buhan (Université Paris Descartes) and Chiara Nardoni (Université Pierre et Marie Curie). Contributors to this project are warmly welcomed.


# Warping
This repository contains the code morphing for warping a source mesh by deforming a closed template mesh. 

### Installation (in progress)
You can grab the sources by cloning this repository or downloading a .zip archive of the sources. In order to build the project, navigate to the created directory and in a command prompt type:

mkdir build

cd build

cmake ..

make

(sudo) make install

### Usage
In a terminal, run:

usage: warping [- nit n] [-t template_file[.mesh]] -s source_file[.mesh]

If the -t template_file[.mesh] is not provided, the code takes in input a default spherical template ("sphere.mesh in the folder demo). 

### Authors and contributors

warping has been initiated by Maya de Buhan (Université Paris Descartes) and Chiara Nardoni (Université Pierre et Marie Curie). Contributors to this project are warmly welcomed.


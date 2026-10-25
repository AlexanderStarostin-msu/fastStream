#Start with the vanila blank Codebade VM

#python wheel for vtk
pip install pyvista vtk

sudo apt update

#headless gnuplot for plotting
sudo apt install gnuplot-nox

#for QGIS
sudo apt install -y gdal-bin libgdal-dev sqlite3 libspatialite-dev python3-pip
pip install geopandas shapely pyogrio trimesh numpy scipy fiona
pip install mapbox-earcut # triangulation for pyvista

#for image generation
sudo apt-get install -y libosmesa6 libosmesa6-dev

#Flex
sudo apt install flex bison

#OpenMPI
sudo apt install -y openmpi-bin libopenmpi-dev
        
#scotch
sudo apt install libscotch-dev
sudo apt install libptscotch-dev

#Ready to build OpenFoam now
cd /workspaces/fastStream/OpenFOAM-10
source etc/bashrc
./Allwmake
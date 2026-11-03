#Start with the vanila blank Codebade VM

logDpd=log.dependenciesInstall
echo "Before installing dependencies" > $logDpd
df -h >> $logDpd

#remove pre-push hook to avoid errors when pushing to remote repo
rm -f .git/hooks/pre-push

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
#simlinks for OpenFOAM-10
sudo ln -sf /usr/lib/x86_64-linux-gnu/libscotch.so /usr/lib/x86_64-linux-gnu/scotch-int32/libscotch-6.1.so
sudo ln -sf /usr/lib/x86_64-linux-gnu/libptscotch.so /usr/lib/x86_64-linux-gnu/scotch-int32/libptscotch-6.1.so
sudo ln -sf /usr/lib/x86_64-linux-gnu/libptscotcherrexit.so /usr/lib/x86_64-linux-gnu/libptscotcherrexit-6.1.so

#Ready to build OpenFoam now
cd /workspaces/fastStream/OpenFOAM-10
source etc/bashrc
./Allwmake

#free space, remove unnecessary packages
sudo apt remove llvm-18
sudo apt remove libspatialite-dev sqlite3
sudo apt autoremove -y

sudo rm -rf /usr/local/sdkman
sudo rm -rf /usr/share/dotnet

sudo apt-get update
sudo apt-get install libboost-dev

echo "After installing dependencies" >> $logDpd
df -h >> $logDpd

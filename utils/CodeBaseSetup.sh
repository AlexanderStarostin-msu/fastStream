#Start with the vanila blank Codebade VM

#python wheels
pip install pyvista vtk

#for image generation
sudo apt-get update
sudo apt-get install -y libosmesa6 libosmesa6-dev

#Flex
sudo apt update
sudo apt install flex bison

#OpenMPI
sudo apt update
sudo apt install -y openmpi-bin libopenmpi-dev
        
#scotch
sudo apt update
sudo apt install libscotch-dev
sudo apt install libptscotch-dev

#Ready to build OpenFoam now
cd OpenFOAM-10
source etc/bashrc
./Allwmake
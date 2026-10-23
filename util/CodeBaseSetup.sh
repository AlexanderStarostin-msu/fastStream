#Start with the vanila blank Codebade VM

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
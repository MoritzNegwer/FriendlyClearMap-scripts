
sudo apt-get update -y
sudo apt-get install -y ca-certificates curl gnupg lsb-release unzip bzip2

#download + unpack Ilastik 
mkdir ./Ilastik
wget -O ./Ilastik/Ilastik.tar.bz2 https://files.ilastik.org/ilastik-1.4.0rc2-Linux.tar.bz2
tar -xf ./Ilastik/Ilastik.tar.bz2

#download + unpack FIJI 
mkdir ./Fiji
wget -O ./Fiji/Fiji.zip https://downloads.imagej.net/fiji/latest/fiji-linux64.zip
unzip ./Fiji/Fiji.zip

#download friendly_clearmap docker image (if available, for testing purposes this is delivered with the script) 
#wget -O ./friendly_clearmap_03.tar ___ 
#assuming it's already here, install it 
sudo docker image load -i ./friendly_clearmap_03.tar 
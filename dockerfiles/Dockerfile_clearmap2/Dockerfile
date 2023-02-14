#import base image
FROM ubuntu:22.04

MAINTAINER Moritz Negwer <moritz.j.negwer@gmail.com>>

#install updates 
RUN apt-get update && apt-get dselect-upgrade -y 

#install python 3.8
RUN apt-get install python3 -y

#install elastix
RUN apt-get install elastix -y

#install wget
RUN apt-get install wget -y

#install libgl.so.1
RUN apt-get install ffmpeg libsm6 libxext6 bzip2 -y

#move to CloudMap folder  
WORKDIR /CloudMap/

#download and unzip Ilastik 
RUN wget --no-check-certificate https://files.ilastik.org/ilastik-1.4.0b13-Linux.tar.bz2
RUN tar xjf ilastik-1.*-Linux.tar.bz2

#download and install anaconda 
RUN wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh 
#install anaconda
RUN bash ./Anaconda* -b

#install gcc
RUN apt-get update 
RUN apt-get dselect-upgrade -y
RUN apt-get install apt-utils -y
RUN apt-get install gcc -y
#RUN apt-get install g++ -y
RUN apt-get install mesa-common-dev -y
RUN apt-get install libgl1-mesa-dev -y
RUN apt-get install libx11-6 -y
RUN apt-get install -y libxi-dev libxcomposite-dev libxfixes-dev libxcursor-dev libxdamage-dev libxrandr-dev libxinerama-dev
RUN apt-get install build-essential -y

#update all anaconda packages
ENV PATH /root/anaconda3/bin:$PATH
RUN echo "y" | conda update --all

#download CM2 yml file (causes docker to crash due to unresolved gtk dependencies, thus manually rebuilt below)
#RUN wget https://christophkirst.github.io/ClearMap2Documentation/html/_downloads/485bb8f940eaf729611aabd5944a1d83/ClearMap.yml


RUN conda install -c conda-forge gtk3

#import modified CM2 yml file
COPY docker_prep/ClearMap.yml ./ClearMap.yml
RUN conda env create -f ./ClearMap.yml

#conda activate environment does not work directly, so we'll do the following instead 
#source: https://pythonspeed.com/articles/activate-conda-dockerfile/

#update all
SHELL ["conda", "run", "-n", "cm2", "/bin/bash", "-c"]
RUN conda update --all

#install h5py
SHELL ["conda", "run", "-n", "cm2", "/bin/bash", "-c"]
RUN conda install -c conda-forge h5py -y

#install pandas
SHELL ["conda", "run", "-n", "cm2", "/bin/bash", "-c"]
RUN conda install pandas -y

#copy in ClearMap2 folder
COPY ClearMap2 ./ClearMap2

#make a symlink for libnifti (broken by recent elastix version)
RUN ln -s /usr/lib/libniftiio.so.0.0.0.0 /usr/lib/libniftiio.so.2

#copy the modified Settings.py over to the FriendlyClearmap folder
COPY docker_prep/Settings.py ./ClearMap2/ClearMap/Settings.py

#copy the modified Plot3d.py to the FriendlyClearmap folder
COPY docker_prep/Plot3d.py ./ClearMap2/ClearMap/Visualization/Plot3d.py

# copy over precompile script
COPY docker_prep/precompile_clearmap2.py ./ClearMap2/ClearMap/precompile_clearmap2.py

#pre-compile clearmap modules 
SHELL ["conda", "run", "-n", "cm2", "/bin/bash", "-c"]
RUN python ./ClearMap2/ClearMap/precompile_clearmap2.py

#fix multiple instances of openmp 
SHELL ["conda", "run", "-n", "cm2", "/bin/bash", "-c"]
RUN conda remove intel-openmp -y
RUN conda install -c conda-forge pytorch -y

#set workdir back to root directory
WORKDIR /

#first command to run once the machine has started up 
#cmd ls /CloudMap/Data/

#run actual CM file 
ENTRYPOINT ["conda", "run", "-n", "cm2", "python", "/CloudMap/Data/process_Ilastik_run_this.py"]


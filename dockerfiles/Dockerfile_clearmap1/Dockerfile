#import base image
FROM ubuntu:20.10

MAINTAINER Moritz Negwer <moritz.j.negwer@gmail.com>>

#install updates 
RUN apt-get update 
RUN apt-get dselect-upgrade -y 

#install python 3.8
RUN apt-get install python -y

#install elastix
RUN apt-get install elastix -y

#install wget
RUN apt-get install wget -y

#install libgl.so.1
RUN apt-get install 'ffmpeg'\
    'libsm6'\ 
    'libxext6'  -y

#move to CloudMap folder  
WORKDIR /CloudMap/

#download and unzip Ilastik 
RUN wget --no-check-certificate https://files.ilastik.org/ilastik-1.4.0b8-Linux.tar.bz2
RUN tar xjf ilastik-1.*-Linux.tar.bz2

#download and install anaconda 
RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh 
#install anaconda
RUN bash ./Anaconda* -b

#update all anaconda packages
ENV PATH /root/anaconda3/bin:$PATH
RUN echo "y" | conda update --all

#download and install clearmap-related things 
RUN conda create --name=cm1 

#conda activate environment does not work directly, so we'll do the following instead 
#source: https://pythonspeed.com/articles/activate-conda-dockerfile/

#SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
#RUN conda install spyder

#install matplotlib
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install matplotlib

#install natsort
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install natsort scipy numpy 

#make following commands run in cm1 environment 
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install tifffile

SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install -c conda-forge imagecodecs-lite

SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install imagecodecs

#install opencv
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install -c conda-forge opencv

#update all
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda update --all

#install pyximport
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install cython

#install gcc
RUN apt-get install gcc -y

#install scikit-image
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install scikit-image

#install h5py
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda install h5py

#update all
SHELL ["conda", "run", "-n", "cm1", "/bin/bash", "-c"]
RUN conda update --all

#copy in FriendlyClearmap folder
COPY FriendlyClearmap ./FriendlyClearmap

#Make Data directory 
#RUN mkdir ./Data/Source/autofluo
#RUN mkdir ./Data/Source/protein

#import data 
#COPY docker_prep/Data ./Data
#(obsolete, will be added as volume when running)

#Make classifiers directory 
RUN mkdir ./classifiers/

#import classifiers (obsolete, will be delivered with data)
#COPY docker_prep/classifiers ./classifiers

#import clearmap scripts (obsolete, will be delivered with data)
#COPY docker_prep/Scripts/ ./FriendlyClearmap/clearmap/ClearMap/Scripts/
 
#make a symlink for libnifti (broken by recent elastix version)
RUN ln -s /usr/lib/libniftiio.so.0.0.0.0 /usr/lib/libniftiio.so.2

#copy the modified Settings.py over to the FriendlyClearmap folder
COPY docker_prep/Settings.py ./FriendlyClearmap/clearmap/ClearMap/Settings.py

#set workdir back to root directory
WORKDIR /

#first command to run once the machine has started up 
#cmd ls /CloudMap/FriendlyClearmap/clearmap/Clearmap/Scripts/docker/

#run actual CM file 
#cmd conda activate cm1 && python /CloudMap/FriendlyClearmap/clearmap/Clearmap/Scripts/docker/process_Ilastik_template.py 
#ENTRYPOINT ["conda", "run", "-n", "cm1", "python", "/CloudMap/FriendlyClearmap/clearmap/Clearmap/Scripts/docker/process_Ilastik_template.py"]
ENTRYPOINT ["conda", "run", "-n", "cm1", "python", "/CloudMap/Data/process_Ilastik_run_this.py"]


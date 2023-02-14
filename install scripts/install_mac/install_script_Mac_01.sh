#! /bin/zsh

#install docker, adapted from here: https://docs.docker.com/desktop/install/mac-install/
#NOTE this assumes you don't have docker installed already. If it's already there, it might get overwritten! 
curl https://desktop.docker.com/mac/main/amd64/Docker.dmg --output ./Docker.dmg

sudo hdiutil attach Docker.dmg
sudo /Volumes/Docker/Docker.app/Contents/MacOS/install
sudo hdiutil detach /Volumes/Docker

#After this, you'll need to confirm the user rights check popup in the graphical user interface. To our knowledge, this cannot be done on the command line via SSH. 



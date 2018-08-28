#!/bin/sh

sudo add-apt-repository main
sudo add-apt-repository universe
sudo add-apt-repository restricted
sudo add-apt-repository multiverse

sudo apt-get update
sudo apt-get autoremove

# enable data mounting
# sudo apt-get -y install open-vm-tools
# sudo apt-get -y install build-essential module-assistant linux-headers-virtual linux-image-virtual && dpkg-reconfigure open-vm-tools
# https://askubuntu.com/questions/29284/how-do-i-mount-shared-folders-in-ubuntu-using-vmware-tools/41386#1051620
sudo mkdir -p /mnt/hgfs
sudo vmhgfs-fuse .host:/ /mnt/hgfs/ -o allow_other -o uid=1000

### AUTO MOUNT in 18.04 server
# Use shared folders between VMWare guest and host
#+goes in /etc/fstab
sudo sh -c "echo '.host:/    /mnt/hgfs/    fuse.vmhgfs-fuse    defaults,allow_other,uid=1000     0    0' >> /etc/fstab"

# install wm tools
# sudo mkdir /mnt/cdrom
# sudo mount /dev/cdrom /mnt/cdrom
# https://kb.vmware.com/s/article/1022525
# tar xzvf /mnt/cdrom/VMwareTools-x.x.x-xxxx.tar.gz -C /tmp/
# #Note: x.x.x-xxxx is the version discovered in the previous step.
# cd /tmp/vmware-tools-distrib/
# sudo ./vmware-install.pl -d

sudo mkdir -p ~/public/phd
export PHD='~/public/phd'
export DATA='/mnt/hgfs/data'
#ln -s ${PHD}/data /mnt/hgfs/data

# install astropy
#sudo apt-get update
#sudo apt-get -y install python-astropy python-numpy python-scipy python-skimage python-photutils

# upgrade photutils to 0.4
#sudo apt-get -y install python-pip
#pip install astropy==2.0.5 --no-deps
#pip install --no-cache-dir --upgrade --no-deps photutils

# install python 3
#sudo apt-get -y install software-properties-common
#sudo add-apt-repository ppa:jonathonf/python-3.6
#sudo apt-get update
#sudo apt-get -y install python3.6

# install python 3 packages 
#sudo apt-get -y install python3-astropy python3-numpy python3-scipy python3-skimage python3-photutils
sudo apt-get update
sudo apt-get -y install python3-pip
python3 -m pip install astropy --no-deps
python3 -m pip install --user numpy scipy matplotlib
python3 -m pip install --no-deps photutils

python3 -m pip install pytest
sudo apt-get -y install python3-tk
python3 -m pip install image_registration

# install iraf
#wget ftp://iraf.noao.edu/iraf/v216/iraf.lnux.x86_64.tar.gzc
#tar zxf iraf.lnux.x86_64.tar.gz
#./install

# install astrometry.net
#sudo apt-get install astrometry.net
#cd /usr/local/astrometry/data
#sudo wget http://data.astrometry.net/4100/index-4110.fits
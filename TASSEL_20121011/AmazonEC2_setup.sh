#!/bin/bash
# A script to configure an EC2 instance for FastqPairedEndToTagCountPlugin analysis
# The script installs the GBS software and Oracle Java 7 (a dependency)

# Use ami-a6ff46cf (CloudBioLinux Ubuntu 12.04 20121019) with 68 Gb RAM (m2.4xlarge),
# and /mnt disk volume to have access to enough disk space for real datasets.
# Run script from /home/ubuntu directory

### NOTE: script assumes that read1.key and read2.key files are packaged with this
###       script in a tgz archive that is uploaded (or emailed) to the EC2 instance
###       and unpacked in the /home/ubuntu directory of the EC2 instance.

# Download data from service lab website directly to EC2 instance, 
# save to /mnt/NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/outputs/fastq
# after this script is run.

sudo chown ubuntu /mnt
cd /mnt
git clone git://github.com/conceptstailored/NCSU-GBS-Pipeline.git
mkdir NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/outputs
mkdir NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/outputs/fastq
mkdir NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/outputs/PairedTBT
cd NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/outputs
### files read1.key and read2.key are in standard TASSEL-GBS format, eg  ###
#Flowcell	Lane	Barcode	Sample	PlateName	Row	Column	Blank
#C116KACXX	7	ACTCCACG	gbs	GBS1	A	1	Blank

# Key files are in archive with this script; create symbolic links from outputs directory to those files

ln -s ~/read1.key read1.key
ln -s ~/read2.key read2.key
chmod 744 /mnt/NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/run_gbs_pairedEnd_pipeline.pl


# Install Java 7 from Oracle, using package from webupd8 ppa;
# see http://www.webupd8.org/2012/01/install-oracle-java-jdk-7-in-ubuntu-via.html
sudo add-apt-repository ppa:webupd8team/java
# Requires a hard return from user to move through choices


sudo apt-get update
sudo apt-get install oracle-java7-installer
# Also requires user input to accept Oracle license terms for Java7, or
sudo echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | sudo /usr/bin/debconf-set-selections

# Install cd-hit clustering software from Debian repository
sudo apt-get install cd-hit

# Run script from outputs directory, giving full path to script after making it executable:
# /mnt/NCSU-GBS-Pipeline/TASSEL_20121011/maizegenetics/run_gbs_pairedEnd_pipeline.pl -Xmx50g -fork1 -FastqPairedEndToTagCountPlugin -i fastq -k read1.key:read2.key -e PstI-MspI:MspI-PstI -s 750000000 -c 10 -o PairedTBT -endPlugin -runfork1

exit


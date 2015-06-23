#!/bin/bash

##### Checking the prerequirement
echo "Checking the system ......"

# check program - blat
BLAT=`type -p blat`
if [ -n "$BLAT" ]; then
	# copy blat to bin directory
	cp $BLAT ./bin
else
	echo "Program 'blat' does not exist! Please download and install 'blatSrc.zip' from <https://users.soe.ucsc.edu/~kent/src/>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# check program - isPcr
ISPCR=`type -p isPcr`
if [ -n "$ISPCR" ]; then
	# copy isPcr to bin directory
	cp $ISPCR ./bin
else
	echo "Program 'isPcr' does not exist! Please download and install 'isPcr.zip' from <https://users.soe.ucsc.edu/~kent/src/>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# check program - primer3_core
PRIMER3=`type -p primer3_core`
if [ -n "$PRIMER3" ]; then
	# copy primer3_core to bin directory
	cp $PRIMER3 ./bin
	# get the directory from primer3_core
	# and copy other programs to ./bin
	PRIMER3_DIR=`dirname $PRIMER3`
	cp "$PRIMER3_DIR/long_seq_tm_test" ./bin
	cp "$PRIMER3_DIR/oligotm" ./bin
	cp "$PRIMER3_DIR/ntdpal" ./bin
	cp "$PRIMER3_DIR/ntthal" ./bin
else
	echo "Program 'primer3_core' does not exist! Please download and install 'Primer3 Release 2.3.6' from <http://primer3.sourceforge.net/releases.php>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# check program - meltsim
MELTSIM=`type -p meltsim`
if [ -n "$MELTSIM" ]; then
	# copy blat to bin directory
	cp $MELTSIM ./bin
else
	echo "Program 'meltsim' does not exist! Please download and install 'MeltSim 1.99.0' from <http://www.bioinformatics.org/meltsim/wiki/>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# check program - cmultiplx
MULTIPLX=`type -p cmultiplx`
if [ -n "$MULTIPLX" ]; then
	# copy cmultiplx to bin directory
	cp $MULTIPLX ./bin
else
	echo "Program 'cmultiplx' does not exist! Please download and install 'MultiPLX version 2.0' from <http://bioinfo.ut.ee/?page_id=167>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# check program - glistmaker
GMASKER=`type -p gmasker`
if [ -n "$GMASKER" ]; then
	# copy primer3_core to bin directory
	cp $GMASKER ./bin
	# get the directory from primer3_core
	# and copy other programs to ./bin
	GMASKER_DIR=`dirname $GMASKER`
	cp "$GMASKER_DIR/glistmaker" ./bin
else
	echo "Program 'gmasker' does not exist! Please download and install 'GenomeMasker package version 1.3' from <http://bioinfo.ut.ee/?page_id=167>."
	echo "Then run the install.sh again. Exiting..."
	exit
fi

# System checking is done
echo "System checking is done!"


##### build the software
echo "Building the software program ......"

# copy bin
cp -r bin ../build
# copy config
cp -r config ../build
# copy examples
cp -r examples ../build
# copy genome
cp -r genome ../build
# copy modules
cp -r modules ../build

# copy mmGTP.conf
cp mmGTP.conf ../build
# copy mouseGTP
cp mouseGTP ../build
# copy mouseGT_PrimerDesign.pl
cp mouseGT_PrimerDesign.pl ../build
# copy mouseGT_SeqPrep.pl
cp mouseGT_SeqPrep.pl ../build

# Software building is done
echo "Program building is done!"


##### install software

# install program in a directory defined by user
read -p "Which directory do you wish to install this program?	" INSTALL_DIR
echo "Installing the program ......" 
cp -r ../build $INSTALL_DIR
mv $INSTALL_DIR/build $INSTALL_DIR/mgtPrimer

# change file permission
chmod 711 $INSTALL_DIR/mgtPrimer/bin/*
chmod 755 $INSTALL_DIR/mgtPrimer/mouseGT*

# create subdirectory
mkdir $INSTALL_DIR/mgtPrimer/tmp
mkdir $INSTALL_DIR/mgtPrimer/result

echo "Program installation is done!"



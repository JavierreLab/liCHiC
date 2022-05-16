#!/bin/bash

# usage
if [ "$#" -ne 2 ]; then
	echo ""
	echo "Usage: 0_mkdir.sh <Work.Dir> <Proj.Name>"
	echo ""
	exit 1
fi

DIR=$1
PROJ=$2

cd $DIR

mkdir $PROJ

cd $PROJ

mkdir data

mkdir data/summary_statistics/
mkdir data/summary_statistics/formatted/
mkdir data/summary_statistics/original/

mkdir jobs/

mkdir logs/

mkdir results/
mkdir results/posterior_probability/
mkdir results/tissue_set_enrichment/


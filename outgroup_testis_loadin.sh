#!/usr/bin/bash

cd /ru-auth/local/home/ulee/scratch_store/analysis/txn_testis_outgroup

declare -a arr=("SRR1024003"
"SRR1024002"
"SRR1617567"
"SRR486109"
"SRR5839446"
"SRR16213299"
"SRR16213300"
"SRR6667440"
"SRR6667439"
"SRR931549"
)


for i in "${arr[@]}"
do
	echo "retrieving $i"
	fasterq-dump $i
done

 # mkdir -p ~/data/scripts ~/data/bin
 # chmod 755 ~/data/scripts ~/data/bin

 # rsync -a rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ~/data/bin/
 # git archive --remote=git://genome-source.soe.ucsc.edu/kent.git \
  # --prefix=kent/ HEAD src/hg/utils/automation \
     # | tar vxf - -C ~/data/scripts --strip-components=5 \
        # --exclude='kent/src/hg/utils/automation/incidentDb' \
      # --exclude='kent/src/hg/utils/automation/configFiles' \
      # --exclude='kent/src/hg/utils/automation/ensGene' \
      # --exclude='kent/src/hg/utils/automation/genbank' \
      # --exclude='kent/src/hg/utils/automation/lastz_D' \
      # --exclude='kent/src/hg/utils/automation/openStack'
  # wget -O ~/data/bin/bedSingleCover.pl 'http://genome-source.soe.ucsc.edu/gitlist/kent.git/raw/master/src/utils/bedSingleCover.pl'
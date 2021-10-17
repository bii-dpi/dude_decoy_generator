#!/bin/csh

set base=$0:h
if ( "$base" == "$0" ) then
    set base="."
endif
set dud = "$base"

awk '{if ($2 == "") print $1,NR,NR; else print $1,NR,$2}' $1 >! input.ism
corina -d neu,rs,wh,rc,mc=1,canon,stergen,preserve -i t=smiles -o t=sdf < input.ism >! corina.sdf
sdconvert -isd corina.sdf -omae corina.mae
@ ok = `epiksqlp.pl lock`
while ($ok == 0)
    epiksqlp.pl
    echo Too many epik licenses checked out... Waiting!
    sleep 10
    @ ok = `epiksqlp.pl lock`
    echo $ok
end
epik -NO_JOBCONTROL -LOCAL -ph 7.0 -pht 1.0 -tp 0.20 -imae corina.mae -omae epik.mae -WAIT
epiksqlp.pl release
run fix_thiazole_mesos.py epik.mae fix.mae
sdconvert -imae fix.mae -osd fix.sdf
convert.py --i fix.sdf --o fix.ism
sed 's/ /\t/g' fix.ism >! tabbed.ism
/raid1/xtalcopy/linux/j2sdk1.4.2_13/jre/bin/java -jar /raid1/soft/mitools/mib.jar -singlepart -prop -onlyOrganic -nlogp -f tabbed.ism -out smi >! mitools.ism
unalias java
set path = ( /usr/arch/jre1.6/bin/ $path )
cxcalc mitools.ism formalcharge | awk 'NR > 1 {print $0}' >! charges.txt
awk '{str1=$1; getline < "charges.txt"; str2=$2; getline < "mitools.ism"; print str1 "\t" $2 "\t" $3 "\t" $7 "\t" $4 "\t" $11 "\t" $8 "\t" $9 "\t" str2}' tabbed.ism >! charged.ism
rm -f corina.sdf corina.mae corina.trc epik.mae epik.log fix.mae fix.sdf fix.ism tabbed.ism mitools.ism charges.txt 
$dud/charged2ligands.py
$dud/decoys.py -s -r -i ligands.charge -o decoys
$dud/translate_user_ids.py
rm -f input.ism charged.ism

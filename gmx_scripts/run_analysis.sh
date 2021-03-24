#!/bin/bash
calc_d='python start_d.py'
gmx=/home/software/science/gromacs-2019_beta/bin/gmx
trajectory=traj.trr
topology=topo.tpr

# Running in each of the replica directories

#molecules
$gmx select -f start.gro  -s $topology    -select "mol 1" "mol 2"  -on  mols.ndx
#molecules backbone l-rmsd
$gmx select -s $topology -on index.ndx -select " mol 1  and name N C CA O" "mol 2 and name N C CA O"
#interface irmsd
$gmx select -f $trajectory  -s  $topology    -select  "((mol 1   and within 1 of mol 2 ) and name N C CA O ) or ((mol 2  and within 1 of mol 1) and name N C CA O)" -on nat_contact.ndx -oi nat_contact.dat -b 0 -e 0
#l-rmsd
echo 0 1 | $gmx rms -f $trajectory -s start.gro -n index.ndx -tu ns -o rms_l_start.xvg -e 100
#i-rmsd
$gmx rms  -s start.gro -f  $trajectory -o irmsd_start.xvg -tu ns -n nat_contact.ndx -e 100
#distance
$gmx distance -s $topology -f $trajectory -select  'com of mol 1 plus com of  mol 2'   -oall -tu ns -pbc yes -e 100
#bsa
echo 1 | $gmx sasa -s $topology -f  $trajectory -o area_all.xvg -e 100 -tu ns
echo 0 | $gmx sasa -s $topology -f $trajectory -o area_1.xvg  -n mols.ndx -e 100 -tu ns
echo 1 | $gmx sasa -s $topology -f  $trajectory -o area_2.xvg  -n mols.ndx -e 100 -tu ns
paste area_all.xvg area_1.xvg area_2.xvg | awk '{print $1, ( $4 + $6 -$2)}' > bsa.xvg
#fnat
echo 0 1 | $gmx hbond -s $topology -f start.gro -n mols.ndx -g hbond_ref.log -hbn hbindex_ref.ndx -contact -r2 0.5
for time in {0..100000..500}
do
echo 0 1 | $gmx hbond -s $topology -f $trajectory  -n mols.ndx -g hbond_"$time".log -hbn hbindex_"$time".ndx -contact -r2 0.5  -e $time -b $time -num num_"$time".xvg
done
#hbond
echo 0 1 | $gmx hbond -s $topology -f $trajectory  -n mols.ndx -tu ns -e 100
$calc_d  -e dist.xvg  -o ddist.xvg
$calc_d  -e bsa.xvg  -o dbsa.xvg
$calc_d  -e fnat_start.txt  -o dfnat_start.xvg

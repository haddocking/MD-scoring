#!/bin/bash
# Running in each of the replica directories

#Change mdps energy groups to have separate ones for proteins and water

sed 's/energygrps               =/energygrps               = mol_1 mol_2  Water /g'  mdout.mdp > mdout.mdp

# Create index files specifying each of the proteins
echo "q" | $gmx make_ndx -f nat_"$x".gro -o index.ndx
cat index.ndx ../mols.ndx > index_mols.ndx
# Rerun trajectories with specified energy groups
gmx grompp -f mdout.mdp -c nat_"$x".gro -n index_mols.ndx -o rer.tpr -p nat_"$x".top -maxwarn 10
gmx mdrun -rerun traj_wat.trr -s rer.tpr

#Extract intermolecular energies from trajectories
echo 22, 23  0 | $gmx energy -f ener.edr -s rer.tpr -o pp.xvg
echo 26, 27, 38, 39 0 | $gmx energy -f ener.edr -s rer.tpr -o  p_wat.xvg
awk '{ print $1 , $2 + $3 }'  pp.xvg > pp_added.xvg
awk '{ print $1 , $2 + $3 + $4 + $5 }'  p_wat.xvg > p_wat_added.xvg


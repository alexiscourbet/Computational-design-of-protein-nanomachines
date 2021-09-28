###acourbet: This script taks an axle and ring, and generate docked configurations by sample rotation and translation along Z

#!/usr/bin/python
import os
import numpy as np
import pyrosetta
import imp
from math import *
from rif.legacy.xyzMath import *
import repeat_utils_w_contact_list 
import hash_subclass_universal


def dock_rotor(p,up_or_down,dist_bins,rot_bins):
    for direction in up_or_down:
      p_up_or_down=repeat_utils_w_contact_list.rotate_around_x(Vec(1.,0.,0.),direction,p.clone())
      for dist in dist_bins:
          p_trans=repeat_utils_w_contact_list.translate_along_z(dist,p_up_or_down.clone())
          for deg in rot_bins:
              print ("Sampling combination",direction, dist,deg)
              c=repeat_utils_w_contact_list.rotate_around_z(Vec(0.,0.,1.0),deg,p_trans.clone())
              yield(c,direction,dist,deg)


#############################################

#init Rosetta
pyrosetta.init()

pose_axle_path="/home/acourbet/rotors/C8D8axle.pdb"
pose_axle=pyrosetta.pose_from_pdb(pose_axle_path)
pose_rotor=pyrosetta.pose_from_pdb("/home/acourbet/rotors/C8D8ring.pdb")

print ("Sampling according to dock params")
dock_gen=dock_rotor(pose_rotor,np.arange(0.,360,360),np.arange(-1,0,1),np.arange(0.0,360,.5))

print ("Outputting pdbs...")

for ipose,direct,distance,degree in dock_gen:
    combined_rotor_clone=pose_axle.clone()
    combined_rotor_clone.append_pose_by_jump(ipose.clone(),1)
    combined_rotor_clone.dump_pdb('sampled_%s_%s_%s.pdb'%(direct,distance,degree))

print ("Done")

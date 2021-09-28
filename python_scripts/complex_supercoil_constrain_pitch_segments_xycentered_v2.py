#!/usr/bin/env python
from os import system,popen
import string
from sys import argv
import math
from math import sin,cos
from xyzMath import *
from scipy.integrate import quad

def enumerate_combinations(list):
    num = len(list)
    ncomb=1
    counter=[]
    for i in range(num):
       ncomb=ncomb*list[i]
       counter.append(0)
    add = 1

    combs=[]
    for j in range(ncomb):
        
      combs.append(string.join(map(lambda x: str(counter[x]), range(num))))  
      for i in range(num):  
        counter[i]=counter[i]+add
        if counter[i]==list[i]:
           counter[i]=0
        else:
           break

    return(combs)

def Generate_Pdbs(Nresl,num_to_output,wl_1,wl_2,Rl_1,Rl_2,Nresl_1,num_link_res,orientation,P0,P1,delta_z,output_file_name,chain_name, chain_order):
 
 deg_to_rad= math.pi/180.
# chain parameters
 
 chain_set=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

 chain_num=chain_set[0:num_to_output]
 
 R1=2.26
#w1 = 720./7.  # 102.85  for a 2-layer heptad repeat
# w1 = 1080./11.
 w1 = 1800./18.
 d=1.51
 z1=0.
 z2=0.
 line='ATOM      7  CA  GLY A   2       5.520   2.352   1.361  1.00 20.00     '
 last=line[54:-1]

 atom_num=1
 res_num=1
 Res_id=[]
 CA_list=[]
 for iter in range(num_to_output):
#  iter=int(chain_order[counter])   
  CA_list.append([])
  Res_id.append([])
#  chain=chain_num[iter]
  chain = chain_name[iter]
  orient = orientation[iter]
  Nres=Nresl[iter]
  Nres_1=Nresl_1[iter]
  if orient == 1:
     res_num=0
  else:
     res_num=Nres+1
  for t in range(Nres+2):      ## need two extra residues to guide placement of 1st and last residue 

#############################################################################
############# First segment #################################################
#############################################################################

    if t <= Nres_1:
     w0_1=wl_1[iter]
     R0_1=Rl_1[iter]
     if orient==1: 
      a1=(w1*(t-1)+P1[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius
     else:
      a1=(w1*t-w1*Nresl[iter]-P1[iter])*deg_to_rad
     alpha=math.asin(R0_1*w0_1*deg_to_rad/d)

     supercoil_phase=P0[iter]+delta_z[iter]*math.tan(alpha)/(R0_1*deg_to_rad)

     a0=(w0_1*(t-1)+supercoil_phase)*deg_to_rad
     if t==Nres_1:
      a0_seg1_last_resi=(w0_1*(Nres_1-1)+supercoil_phase)*deg_to_rad
     x=R0_1*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
     y=R0_1*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
     if w0_1==0:
      z=d*(t-1)+delta_z[iter]
      z1=d*(t-1)+delta_z[iter]
      if t==Nres_1:
       z_seg1_last_res=d*(Nres_1-1)+delta_z[iter]
     else:
      z=R0_1*w0_1*(t-1)*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+delta_z[iter]
      z1= R0_1*w0_1*(t-1)*deg_to_rad/math.tan(alpha)+delta_z[iter]
      if t==Nres_1:
       z_seg1_last_res=R0_1*w0_1*(Nres_1-1)*deg_to_rad/math.tan(alpha)+delta_z[iter]
     if t==1:
      z1_first_res=delta_z[iter]
#     print t, z1, a0/3.1415926*180
##############################################################################
############# Second segment #################################################
#############################################################################

    elif t >= (Nres_1 + 1 + num_link_res[iter]):
     w0_2=wl_2[iter]
#      w0_2=wl[iter]/math.sqrt((1+(math.pow(Rl_2[iter],2)-math.pow(Rl_1[iter],2))*math.pow(wl[iter]*math.pi/180,2)/math.pow(1.51,2)))
     R0_2=Rl_2[iter]
     if orient==1: 
      a1=(w1*(t-1)+P1[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius
     else:
      a1=(w1*t-w1*Nresl[iter]-P1[iter])*deg_to_rad
     alpha=math.asin(R0_2*w0_2*deg_to_rad/d)

#     a0=(w0_2*(t-1)+supercoil_phase)*deg_to_rad
     a0=(w0_2*(t-(Nres_1+num_link_res[iter]))+a0_link_last_resi/deg_to_rad)*deg_to_rad
#     a0=w0_2*deg_to_rad+a0

     x=R0_2*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
     y=R0_2*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
     if w0_2==0:
      z=d*(t-(Nres_1+num_link_res[iter]))+z_link_last_res
      z1=d*(t-(Nres_1+num_link_res[iter]))+z_link_last_res
      if t==Nres:
       z1_last_res=d*(Nres-(Nres_1+num_link_res[iter]))+z_link_last_res
     else:
      z= R0_2*w0_2*(t-(Nres_1+num_link_res[iter]))*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+z_link_last_res
      z1= R0_2*w0_2*(t-(Nres_1+num_link_res[iter]))*deg_to_rad/math.tan(alpha)+z_link_last_res
      if t==Nres:
       z1_last_res=R0_2*w0_2*(Nres-(Nres_1+num_link_res[iter]))*deg_to_rad/math.tan(alpha)+z_link_last_res
#     print t, z1, a0/3.1415926*180
##############################################################################
############# Linker segment #1 ##############################################
##############################################################################

    elif Nres_1 < t <= Nres_1+int(num_link_res[iter]/2):
     R0_l=Rl_1[iter]+(Rl_2[iter]-Rl_1[iter])/num_link_res[iter]*(t-Nres_1)
#     R0_l_p=Rl_1[iter]+(Rl_2[iter]-Rl_1[iter])/num_link_res[iter]*(t-1-Nres_1)
#     def integrand(R_l):
#      return -math.pow(-3.18/180*math.pi,3)/math.pow(1.51,2)*R_l*math.pow((1+(math.pow(R_l,2)-math.pow(9.5,2))*math.pow(-3.18/180*math.pi,2)/math.pow(1.51,2)),-1.5)
#     delta_w0_l,err= quad(integrand, R0_l_p,R0_l)
#     if t==Nres_1+1:
#      w0_l=w0_1+delta_w0_l*180/math.pi
#     else:
#      w0_l=w0_l+delta_w0_l*180/math.pi
#     w0_l=integrate((-R_l/math.pow(1.51,2)*(math.pow(math.pow(wl[iter],-2)+(math.pow(R_l,2)-math.pow(Rl_1[iter],2))/math.pow(1.51,2),-1.5))),(R_l,R0_l,R0_l_p))
     w0_l=wl_1[iter]/math.sqrt((1+(math.pow(R0_l,2)-math.pow(Rl_1[iter],2))*math.pow(wl_1[iter]*math.pi/180,2)/math.pow(1.51,2)))
     if orient==1: 
      a1=(w1*(t-1)+P1[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius
     else:
      a1=(w1*t-w1*Nresl[iter]-P1[iter])*deg_to_rad
     alpha=math.asin(R0_l*w0_l*deg_to_rad/d)

#     a0=w0_l*deg_to_rad+a0
#     a0=(w0_l*(t-1)+supercoil_phase)*deg_to_rad
     a0=(w0_l*(t-Nres_1)+a0_seg1_last_resi/deg_to_rad)*deg_to_rad
     x=R0_l*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
     y=R0_l*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
     if w0_l==0:
      z=d*(t-Nres_1)+z_seg1_last_res
      z1=d*(t-Nres_1)+z_seg1_last_res
     else:
      z= R0_l*w0_l*(t-Nres_1)*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+z_seg1_last_res
      z1= R0_l*w0_l*(t-Nres_1)*deg_to_rad/math.tan(alpha)+z_seg1_last_res
     if t == Nres_1+int(num_link_res[iter]/2):
      t_turn=t
      z_turn=z1
      a0_turn=a0
#      print "z_turn:",z_turn,"a0_turn:",a0_turn/3.1415926*180
#     print t, z1, a0/3.1415926*180
    
##############################################################################
############# Linker segment #2 ##############################################
##############################################################################
    else:
     R0_l=Rl_1[iter]+(Rl_2[iter]-Rl_1[iter])/num_link_res[iter]*(t-Nres_1)
     w0_l=wl_2[iter]/math.sqrt((1+(math.pow(R0_l,2)-math.pow(Rl_2[iter],2))*math.pow(wl_2[iter]*math.pi/180,2)/math.pow(1.51,2)))
     if orient==1: 
      a1=(w1*(t-1)+P1[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius
     else:
      a1=(w1*t-w1*Nresl[iter]-P1[iter])*deg_to_rad
     alpha=math.asin(R0_l*w0_l*deg_to_rad/d)

     a0=(w0_l*(t-t_turn)+a0_turn/deg_to_rad)*deg_to_rad
     if t==Nres_1+num_link_res[iter]:
      a0_link_last_resi=(w0_l*(Nres_1+num_link_res[iter]-t_turn)+a0_turn/deg_to_rad)*deg_to_rad
     x=R0_l*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
     y=R0_l*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
     if w0_l==0:
      z=d*(t-t_turn)+z_turn
      z1=d*(t-t_turn)+z_turn
      if t==Nres_1+num_link_res[iter]:
       z_link_last_res=d*(Nres_1+num_link_res[iter]-t_turn)+z_turn
     else:
      z= R0_l*w0_l*(t-t_turn)*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+z_turn
      if t==Nres_1+num_link_res[iter]:
       z_link_last_res= R0_l*w0_l*(Nres_1+num_link_res[iter]-t_turn)*deg_to_rad/math.tan(alpha)+z_turn
      z1= R0_l*w0_l*(t-t_turn)*deg_to_rad/math.tan(alpha)+z_turn
#     print t, z1, a0/3.1415926*180 
    


    CA_list[iter].append( (res_num,Vec(x,y,z)) )
    Res_id[iter].append(res_num)
    atom_num=atom_num+1
     
    if orient  == 1:
      res_num=res_num+1
    else:
      res_num=res_num-1

#######this will place the middle of helix to xy plane
  for i in range(len(CA_list[iter])):
    CA_list[iter][i][1].z=CA_list[iter][i][1].z-(z1_last_res-z1_first_res)/2
     

 # convert CA trace to full backbone model by superimposing on ideal template
 # by matching 3 consecutive CA atoms
 # set up ideal template
 stub_file=map(string.split,open('ideal.pdb','r').readlines())
 atom=[]
 for line in stub_file:
    atom.append( (Vec(float(line[6]),float(line[7]),float(line[8]))))
    
 ideal_stub=stub(atom[6],atom[1],atom[11])
##  print atom[6]
##  print atom[8]
##  print Vec.distance(atom[6],atom[8]),'dist'
#now make full backbone pdb
 full_pdb=open('%s'%(output_file_name),'w')
 atom_num=1

 res_num=0 
 for counter in range(num_to_output):
  iter=int(chain_order[counter])
  chain=chain_name[iter]
  CA_chain_u=CA_list[iter]
  CA_chain = sorted(CA_chain_u, key = lambda res: res[0])
  for res in range(1,Nresl[iter]+1):
    res_num=res_num+1
#    print res, res_num
#    print res_num
#    print CA_chain[res][1]
    actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])
    transform=actual_stub * ~ideal_stub

#    start_d=Vec.distance(atom[5],atom[6])
#    end_d = Vec.distance(coords,transform*atom[6])
#    print start_d,end_d,'dist'
#    print CA_list[res],'ori',coords

# N
    coords=transform*atom[5]
    full_pdb.write('ATOM %6d  N   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

# CA   (use actual CA from trace rather than superimposed one)
    coords=CA_chain[res][1]
    tcoords=transform*atom[6]
#    print coords,tcoords,'CA'
    
    full_pdb.write('ATOM %6d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

#  NH
    coords=transform*atom[7]
    full_pdb.write('ATOM %6d  H   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

#  C
    coords=transform*atom[8]
    full_pdb.write('ATOM %6d  C   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

# O
    coords=transform*atom[9]
    full_pdb.write('ATOM %6d  O   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
    atom_num=atom_num+1

    start_d=Vec.distance(atom[8],atom[6])
    end_d = Vec.distance(transform*atom[8],transform*atom[6])
#    print start_d,end_d,'dist C  CA '
    
 return()

def generate_values(data_range):
    
        inc=(data_range[1]-data_range[0])/max(data_range[2]-1.,1.)
        vals=[]          
        for j in range(data_range[2]):
              vals.append(data_range[0]+inc*j)
        return(vals)
    
def input_params(input):

    output_file_prefix=input[0][0]

    Nres=[]
    w0l_1=[]
    w0l_2=[]
    R0l_1=[]
    R0l_2=[]
    Nres_1=[]
    num_link_res=[]
    P0l =[]
    orientation=[]
    P1l = []
    Zl = []
    num_to_output=int(input[1][0])
    chain_name=[]
    chain_order=[]
    for i in range(num_to_output):
        chain_name.append(input[2][i])
        chain_order.append(input[3][i])
    num_vars=11
    for i in range(num_to_output):
        offset=num_vars*(i)+4
        Nres.append(int(input[0+offset][0]))
        Nres_1.append(int(input[1+offset][0]))
        num_link_res.append(int(input[2+offset][0]))
        w0l_1.append( (float(input[3+offset][0]),float(input[3+offset][1]),int(input[3+offset][2])))
        w0l_2.append( (float(input[4+offset][0]),float(input[4+offset][1]),int(input[4+offset][2])))
        R0l_1.append( (float(input[5+offset][0]),float(input[5+offset][1]),int(input[5+offset][2])))
        R0l_2.append( (float(input[6+offset][0]),float(input[6+offset][1]),int(input[6+offset][2])))
        P0l.append( (float(input[7+offset][0]),float(input[7+offset][1]),int(input[7+offset][2])))
        orientation.append( int(input[8+offset][0] ))
        P1l.append( (float(input[9+offset][0]),float(input[9+offset][1]),int(input[9+offset][2])))
        Zl.append( (float(input[10+offset][0]),float(input[10+offset][1]),int(input[10+offset][2])))

        
    print ' ###############    SUPERHELIX PARAMS   ############## \n'
    print ' output file prefix:  %s \n'%(output_file_prefix)
    print ' number of chains:  %s \n'%(num_to_output)
    for i in range(num_to_output):
       print 'HELIX %s PARAMS ############################################### \n'%(i)           
       print ' helix length, orientation, and id:  %s %s %s \n'%(Nres[i],orientation[i],chain_name[i])
       print ' starting twist for 1st segment: %s  ending twist: %s  number of samples: %s \n'%(w0l_1[i][0],w0l_1[i][1],w0l_1[i][2])
       print ' starting twist for 2nd segment:: %s  ending twist: %s  number of samples: %s \n'%(w0l_2[i][0],w0l_2[i][1],w0l_2[i][2])
       print ' starting R0 for 1st segment:   %s   ending R0:  %s   number of samples:  %s \n'%(R0l_1[i][0],R0l_1[i][1],R0l_1[i][2])
       print ' starting R0 for 2nd segment:   %s   ending R0:  %s   number of samples:  %s \n'%(R0l_2[i][0],R0l_2[i][1],R0l_2[i][2])
       print ' starting P0:   %s   ending P0:  %s   number of samples:  %s \n'%(P0l[i][0],P0l[i][1],P0l[i][2])
       print ' starting P1 for 1st segment:   %s   ending P1:  %s   number of samples:  %s \n'%(P1l[i][0],P1l[i][1],P1l[i][2])                 
       print ' starting Z:   %s   ending Z:  %s   number of samples:  %s \n'%(Zl[i][0],Zl[i][1],Zl[i][2])                 

    print ' ##########  HELIX ORIENTATION CHAIN CHAIN_ORDER ###############  \n '
    for i in range(num_to_output):
        print '%s %s %s %s  \n'%(i, orientation[i],chain_name[i],chain_order[i])
    # sample p1 evenly, w and R around the input values

    ## z_list=[1.5,2.0,2.5,3.0]
    ## p1_s=125.
    ## p2_s=235.


    number_of_combinations=1
    for i in range(num_to_output):

 
#        number_of_combinations=number_of_combinations*R0l_1[i][2]*R0l_2[i][2]*P0l[i][2]*P1l[i][2]*P1l_2[i][2]*Zl[i][2]
        if w0l_1[i][2] > 0. : number_of_combinations=number_of_combinations*w0l_1[i][2]
        if w0l_2[i][2] > 0. : number_of_combinations=number_of_combinations*w0l_2[i][2]
        if R0l_1[i][2] > 0. : number_of_combinations=number_of_combinations*R0l_1[i][2]
        if R0l_2[i][2] > 0. : number_of_combinations=number_of_combinations*R0l_2[i][2]
        if P0l[i][2] > 0. : number_of_combinations=number_of_combinations*P0l[i][2]
        if P1l[i][2] > 0. : number_of_combinations=number_of_combinations*P1l[i][2]
        if Zl[i][2] > 0. : number_of_combinations=number_of_combinations*Zl[i][2]

    print  '  NUMBER OF PDBS TO GENERATE:  %s \n'%number_of_combinations

# mofify below
    combinations=[]
    w0_1=[]
    w0_2=[]
    R0_1=[]
    R0_2=[]
    P0=[]
    P1=[]
    Z=[]
    for i in range(num_to_output):

 
        if w0l_1[i][2] > 0.:
            w0_1.append(generate_values(w0l_1[i]))
        else:
            w0_1.append(["constrained"])
        if w0l_2[i][2] > 0.:
            w0_2.append(generate_values(w0l_2[i]))
        else:
            w0_2.append(["constrained"])
        if R0l_1[i][2] > 0.:
            R0_1.append(generate_values(R0l_1[i]))
        else:
            R0_1.append(["constrained"])
        if R0l_2[i][2] > 0.:
            R0_2.append(generate_values(R0l_2[i]))
        else:
            R0_2.append(["constrained"])
        if P1l[i][2] > 0.:
            P1.append(generate_values(P1l[i]))
        else:
            P1.append(["constrained"])
        P0.append(generate_values(P0l[i]))
        Z.append(generate_values(Zl[i]))

             
    return(Nres,num_to_output,orientation,R0_1,R0_2,Nres_1,num_link_res,w0_1,w0_2,P0,P1,Z,output_file_prefix,chain_name,chain_order)

####################################################                      

input_file=argv[1]
tag = argv[2]
input =map(string.split,open(input_file,'r').readlines())
Nres,num_to_output,orientation,R0_1,R0_2,Nres_1,num_link_res,w0_1,w0_2,P0,P1,Z,output_file_prefix,chain_name,chain_order=input_params(input)
items=[]
for i in range(num_to_output):
  items.append(len(w0_1[i]))
  items.append(len(w0_2[i]))
  items.append(len(R0_1[i]))
  items.append(len(R0_2[i]))
  items.append(len(P0[i]))
  items.append(len(P1[i]))
  items.append(len(Z[i]))
combos=enumerate_combinations(items)
for combo in combos:
    id=map(int,string.split(combo))

    w0v_1=[]
    w0v_2=[]
    R0v_1=[]
    R0v_2=[]
    P0v=[]
    P1v=[]
    Zv=[]
    for i in range(num_to_output):
        if len(w0_1[i])==1 and w0_1[i][0]=="constrained" and R0_1[i][0]=="constrained":
            w0v_1.append(w0_1[0][id[0]])
        elif len(w0_1[i])==1 and w0_1[i][0]=="constrained" and R0_1[i][0]!="constrained":
            w0v_1.append( w0_1[0][id[0]]/math.sqrt(1+(math.pow(R0_1[i][id[2+i*7]],2)-math.pow(R0_1[0][id[2]],2))*math.pow(w0_1[0][id[0]]*math.pi/180,2)/math.pow(1.51,2)))
        else:
            w0v_1.append(w0_1[i][id[0+i*7]])
        if len(w0_2[i])==1 and w0_2[i][0]=="constrained" and R0_2[i][0]=="constrained":
            w0v_2.append(w0_2[0][id[1]])
        elif len(w0_2[i])==1 and w0_2[i][0]=="constrained" and R0_2[i][0]!="constrained":
	    w0v_2.append( w0_2[0][id[1]]/math.sqrt(1+(math.pow(R0_2[i][id[3+i*7]],2)-math.pow(R0_1[0][id[3]],2))*math.pow(w0_2[0][id[1]]*math.pi/180,2)/math.pow(1.51,2)))
        else:
            w0v_2.append(w0_2[i][id[1+i*7]])
        if len(R0_1[i])==1 and R0_1[i][0]=="constrained":
            R0v_1.append(R0_1[0][id[2]])
        else:
            R0v_1.append(R0_1[i][id[2+i*7]])
        if len(R0_2[i])==1 and R0_2[i][0]=="constrained":
            R0v_2.append(R0_2[0][id[3]])
        else:
            R0v_2.append(R0_2[i][id[3+i*7]])
        if len(P1[i])==1 and P1[i][0]=="constrained":
            P1v.append(P1[0][id[5]])
        else:
            P1v.append(P1[i][id[5+i*7]])
        P0v.append(P0[i][id[4+i*7]])
        Zv.append(Z[i][id[6+i*7]])           
        
##     ## to get close to C2 symmetry for axis in x-y plane, switch sign of  phase and delta_z of sym mates
##         ## 
## ##     if num_chain==2:
## ## 	helix_phase.append(orientation[1]*helix_phase[0])
## ##         delta_z.append(orientation[1]*delta_z[0])

## ##     if num_chain==4:   ## warning-this is 4 helix bundle centric
## ## 	helix_phase.append(-helix_phase[1])
## ##         helix_phase.append(-helix_phase[0])
## ##         delta_z.append(-delta_z[1])
## ##         delta_z.append(-delta_z[0])

 
    out_file_name='%s'%(tag)
    for i in range(num_to_output):
      out_file_name=out_file_name+'_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f'%(w0v_1[i],w0v_2[i],R0v_1[i],R0v_2[i],P0v[i],P1v[i],Zv[i])
#       out_file_name=out_file_name
    
    out_file_name=out_file_name+'.pdb'
    Generate_Pdbs(Nres,num_to_output,w0v_1,w0v_2,R0v_1,R0v_2,Nres_1,num_link_res,orientation,P0v,P1v,Zv,out_file_name,chain_name,chain_order)
    
## for i in range(3):
##     p1=p1_s+(i-1)*5
##     for ii in range(3):
##      p2=p2_s+(ii-1)*5   
##      for j in range(3):
##         R=R0+(j-1)*0.1
##         for k in range(5):
##             w=w0+(k-2)*0.15
##             for l in range(4):
##                 delta_z=z_list[l]
##                 output_file='fine2_%s_%s_%s_%s_%s_%s'%(w,R,p1,p2,delta_z,output_file_name)
## #    print output_file_name,Nres,w,R,p1
##                 Generate_Pdbs(Nres,w,R,num_chain,p1,p2,delta_z,output_file,num_to_output)
## #print atom_list

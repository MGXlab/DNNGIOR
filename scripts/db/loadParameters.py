# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 12:28:47 2022

@author: danie
"""

from pathlib import Path
import os
import sys

sys.path.append(os.path.join(Path(os.getcwd()).parents[0], 'db'))

from updateParameters import *
from lmfit import Parameters



def getPramsFromFile(species, filePath):
    lmfit_params = Parameters()    
    with open(filePath) as f:
        f.readline()
        for line in f:
            a = line.strip().split('\t')
            if a[0] == species:
                lmfit_params.add(a[1], value = float(a[2]), min = float(a[3]), max = float(a[4]), vary=True)
    
    return lmfit_params
    
    


def assignBhParams(lmfit_params, conn):
    
    a1 = str(lmfit_params['bh_z1_t1'].value**(lmfit_params['bh_z1_h1'].value))
    
    b1 = "metObj.metD['trehalose'].concentration**" + str(lmfit_params['bh_z1_h1'].value)
    
    a2 = str(lmfit_params['bh_z1_t2'].value**(lmfit_params['bh_z1_h2'].value))
    
    b2 = "metObj.metD['glucose'].concentration**" + str(lmfit_params['bh_z1_h2'].value)

        
    #zeta1 = "(" + a1 + "/(" + a1 + " + " + b1 + "))"#"*(" + b2 + "/(" + a2 + " + " + b2 + "))"
    zeta1 = "metObj.metD['trehalose'].concentration < " +  str(lmfit_params['bh_z1_t1'].value)
    
    zeta2 = '""'
    
    a1 = str(lmfit_params['bh_z3_t1'].value**(lmfit_params['bh_z3_h1'].value))
    
    b1 = "metObj.metD['glucose'].concentration**" + str(lmfit_params['bh_z3_h1'].value)
    
    a2 = str(lmfit_params['bh_z3_t2'].value**(lmfit_params['bh_z3_h2'].value))
    
    b2 = "metObj.metD['glutamate'].concentration**" + str(lmfit_params['bh_z3_h2'].value)
    
    
    #zeta3 = "(((" + a1 + ")/(" + a1 + "+ " + b1 + ")) + ((" + a2 + ")/(" + a2 + " + " + b2 + ")) + ((((" + a1 + ")/(" + a1 + " + " + b1 + ")) - ((" + a2 + ")/(" + a2 + " + " + b2 + ")))**2)**0.5)/2"
    
    zeta3 = "metObj.metD['glucose'].concentration < " + str(lmfit_params['bh_z3_t1'].value) + " or metObj.metD['glutamate'].concentration < " + str(lmfit_params['bh_z3_t2'].value)
    
    zeta4 = '""'
    
    with conn:
        
        update_subpopulations(conn, (lmfit_params['bh_expA_mumax'].value, lmfit_params['bh_expA_pHopt'].value, lmfit_params['bh_expA_pHalpha'].value, 'bh.expA'))
        
        update_subpopulations(conn, (lmfit_params['bh_expB_mumax'].value, lmfit_params['bh_expB_pHopt'].value, lmfit_params['bh_expB_pHalpha'].value, 'bh.expB'))
        
        
        
        update_subpopulations2subpopulations(conn, (zeta1, lmfit_params['bh_z1_r'].value, 6))
        
        update_subpopulations2subpopulations(conn, (zeta2, lmfit_params['bh_z2_r'].value, 7))
        
        update_subpopulations2subpopulations(conn, (zeta3, lmfit_params['bh_z3_r'].value, 8))
        
        update_subpopulations2subpopulations(conn, (zeta4, lmfit_params['bh_z4_r'].value, 9))
        
        
        
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expA_ft1_trehalose_g'].value, lmfit_params['bh_monod_trehalose'].value, 24))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expA_ft1_acetate_g'].value, 0, 25))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expA_ft1_lactate_g'].value, 0, 26))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expA_ft2_pyruvate_g'].value, lmfit_params['bh_monod_pyruvate'].value, 27))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expA_ft2_acetate_g'].value, 0, 28))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expB_ft3_glutamate_g'].value, lmfit_params['bh_monod_glutamate'].value, 29))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expB_ft3_glucose_g'].value, lmfit_params['bh_monod_glucose'].value, 30))
        
        update_feedingTerms2metabolites(conn, (lmfit_params['bh_expB_ft3_acetate_g'].value, 0, 31))
        
        update_wc(conn, (lmfit_params['glutamate_conc'].value, "glutamate"))
        
    
    

    
    

def assignParameters2db(species, lmfit_params, conn):
    
    
    
    
    
    if species == 'bh':
        assignBhParams(lmfit_params, conn)
    elif species == 'bt':
        
        zeta1_bt = "((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['bt_z1_h'].value) + ") /((" + str(lmfit_params['bt_z1_t'].value**lmfit_params['bt_z1_h'].value) + ") + ((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['bt_z1_h'].value) + "))"
       
        
        zeta2_bt ="(" + str(lmfit_params['bt_z2_t'].value**lmfit_params['bt_z2_h'].value) + ")/((" + str(lmfit_params['bt_z2_t'].value**lmfit_params['bt_z2_h'].value) + ") + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " + str(lmfit_params['bt_z2_h'].value) + ")"
        
        
        zeta3_bt = "(" + str(lmfit_params['bt_z3_t1'].value**lmfit_params['bt_z3_h1'].value) + ")/((" + str(lmfit_params['bt_z3_t1'].value**lmfit_params['bt_z3_h1'].value) + ") + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " + str(lmfit_params['bt_z3_h1'].value) + ") * (metObj.metD['mannose'].concentration ** " + str(lmfit_params['bt_z3_h2'].value) + ") / ((" + str(lmfit_params['bt_z3_t2'].value**lmfit_params['bt_z3_h2'].value) + ") + (metObj.metD['mannose'].concentration ** " + str(lmfit_params['bt_z3_h2'].value) + "))"
        
        
        zeta4_bt = "(" + str(lmfit_params['bt_z4_t'].value**lmfit_params['bt_z4_h'].value) + ")/((" + str(lmfit_params['bt_z4_t'].value**lmfit_params['bt_z4_h'].value) + ") + (metObj.metD['mannose'].concentration) ** " + str(lmfit_params['bt_z4_h'].value) + ")"
        
        zeta5_bt = '""'
       
        with conn:
            
            

            update_subpopulations(conn, (lmfit_params['bt_lag_mumax'].value, 7.0, 5, 'bt.lag'))
            update_subpopulations(conn, (lmfit_params['bt_expA_mumax'].value, lmfit_params['bt_expA_pHopt'].value, lmfit_params['bt_expA_pHalpha'].value, 'bt.expA'))
            update_subpopulations(conn, (lmfit_params['bt_expB_mumax'].value, lmfit_params['bt_expB_pHopt'].value, lmfit_params['bt_expB_pHalpha'].value, 'bt.expB'))
            
            update_subpopulations2subpopulations(conn, (zeta1_bt, lmfit_params['bt_z1_r'].value, 1))
            update_subpopulations2subpopulations(conn, (zeta2_bt, lmfit_params['bt_z2_r'].value, 2))
            update_subpopulations2subpopulations(conn, (zeta3_bt, lmfit_params['bt_z3_r'].value, 3))
            update_subpopulations2subpopulations(conn, (zeta4_bt, lmfit_params['bt_z4_r'].value, 4))
            update_subpopulations2subpopulations(conn, (zeta5_bt, lmfit_params['bt_z5_r'].value, 5))
            
            
            
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft1_glucose_g'].value, lmfit_params['bt_monod_glucose'].value, 1))
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft1_succinate_g'].value, 0, 2))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft1_acetate_g'].value, 0, 3))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft1_lactate_g'].value, 0, 4))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft1_formate_g'].value, 0, 5))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft2_pyruvate_g'].value, lmfit_params['bt_monod_pyruvate'].value, 6))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft2_acetate_g'].value, 0, 7))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft2_lactate_g'].value, 0, 8))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_lag_ft2_formate_g'].value, 0, 9))
            

            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft1_glucose_g'].value, lmfit_params['bt_monod_glucose'].value, 10))
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft1_succinate_g'].value, 0, 11))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft1_acetate_g'].value, 0, 12))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft1_lactate_g'].value, 0, 13))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft1_formate_g'].value, 0, 14))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft2_pyruvate_g'].value, lmfit_params['bt_monod_pyruvate'].value, 15))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft2_acetate_g'].value, 0, 16))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft2_lactate_g'].value, 0, 17))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expA_ft2_formate_g'].value, 0, 18))
            
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expB_ft3_mannose_g'].value, lmfit_params['bt_monod_mannose'].value, 19))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expB_ft3_succinate_g'].value, 0, 20))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expB_ft3_acetate_g'].value, 0, 21))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expB_ft3_lactate_g'].value, 0, 22))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['bt_expB_ft3_formate_g'].value, 0, 23))
            
            
    elif species == 'ri':
        zeta1_ri = "((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['ri_z1_h'].value) + ") /((" + str(lmfit_params['ri_z1_t'].value**lmfit_params['ri_z1_h'].value) + ") + ((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['ri_z1_h'].value) + "))"
        


        zeta2_ri = "(metObj.metD['lactate'].concentration ** " +  str(lmfit_params['ri_z2_h'].value) + ") /((" + str(lmfit_params['ri_z2_t'].value**lmfit_params['ri_z2_h'].value) + ") + (metObj.metD['lactate'].concentration ** " +  str(lmfit_params['ri_z2_h'].value) + "))"
        
        
     

        zeta3_ri = "(" + str(lmfit_params['ri_z3_t'].value**lmfit_params['ri_z3_h'].value) + ")/((" + str(lmfit_params['ri_z3_t'].value**lmfit_params['ri_z3_h'].value) + ") + (metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " + str(lmfit_params['ri_z3_h'].value) + ")"
        
        
        zeta4_ri = '""'
        
        zeta5_ri = '""'
        
        zeta6_ri = "((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['ri_z6_h'].value) + ") /((" + str(lmfit_params['ri_z6_t'].value**lmfit_params['ri_z6_h'].value) + ") + ((metObj.metD['glucose'].concentration + metObj.metD['pyruvate'].concentration) ** " +  str(lmfit_params['ri_z6_h'].value) + "))"
        
       
        with conn:
            
            update_subpopulations(conn, (lmfit_params['ri_lag_mumax'].value, 7.0, 5.0, 'ri.lag'))
            
            update_subpopulations(conn, (lmfit_params['ri_expA_mumax'].value, lmfit_params['ri_expA_pHopt'].value, lmfit_params['ri_expA_alpha'].value, 'ri.expA'))
             
            update_subpopulations(conn, (lmfit_params['ri_maint_mumax'].value, lmfit_params['ri_maint_pHopt'].value, lmfit_params['ri_maint_alpha'].value, 'ri.maint'))
            
            
            update_subpopulations2subpopulations(conn, (zeta1_ri, lmfit_params['ri_z1_r'].value, 10))
            
            update_subpopulations2subpopulations(conn, (zeta2_ri, lmfit_params['ri_z2_r'].value, 11))
            
            update_subpopulations2subpopulations(conn, (zeta3_ri, lmfit_params['ri_z3_r'].value, 12))
            
            update_subpopulations2subpopulations(conn, (zeta4_ri, lmfit_params['ri_z4_r'].value, 13))
            
            update_subpopulations2subpopulations(conn, (zeta5_ri, lmfit_params['ri_z5_r'].value, 14))
            
            update_subpopulations2subpopulations(conn, (zeta5_ri, lmfit_params['ri_z6_r'].value, 15))
            
            
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft1_glucose_g'].value, lmfit_params['ri_monod_glucose'].value, 37))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft1_lactate_g'].value, 0.0, 38))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft1_acetate_g'].value, 0.0, 39))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft1_butyrate_g'].value, 0.0, 40))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft2_pyruvate_g'].value, lmfit_params['ri_monod_pyruvate'].value, 41))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft2_acetate_g'].value, 0.0, 42))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_lag_ft2_butyrate_g'].value, 0.0, 43))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft1_glucose_g'].value, lmfit_params['ri_monod_glucose'].value, 44))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft1_lactate_g'].value, 0.0, 45))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft1_acetate_g'].value, 0.0, 46))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft1_butyrate_g'].value, 0.0, 47))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft2_pyruvate_g'].value, lmfit_params['ri_monod_pyruvate'].value, 48))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft2_acetate_g'].value, 0.0, 49))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_expA_ft2_butyrate_g'].value, 0.0, 50))
            
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_maint_ft3_lactate_g'].value, lmfit_params['ri_monod_lactate'].value, 51))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_maint_ft3_butyrate_g'].value, 0.0, 52))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_maint_ft4_lactate_g'].value, lmfit_params['ri_monod_lactate'].value, 53))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_maint_ft4_acetate_g'].value, lmfit_params['ri_monod_actetate'].value, 54))
            
            update_feedingTerms2metabolites(conn, (lmfit_params['ri_maint_ft4_butyrate_g'].value, 0.0, 55))
            
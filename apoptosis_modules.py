# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:27:29 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

def apoptosis_monomers():
    # **Activators**.
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})

    # **Effectors**
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', ['bf', 's1', 's2', 'state', 'S184'], {'state':['C', 'M', 'A'], 'S184':['U', 'P']})

    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']})

    # **Anti-Apoptotics**
    Monomer('Bcl2', ['bf', 'state', 'S70'], {'state': ['C', 'M'], 'S70':['U', 'P']})
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', ['bf', 'state'], {'state':['M', 'C']})

    # **Sensitizers**
    Monomer('Bad', ['bf', 'state', 'S75', 'S99', 'S112'], {'state':['C', 'M'], 'S75': ['U', 'P'], 'S99':['U', 'P'], 'S112':['U', 'P']}) #S75: phosphorylation site for S6K; S99: major site of AKT phosphorylation
    Monomer('Noxa', ['bf', 'state'], {'state': ['C', 'M']})
    Monomer('Bim', ['bf', 'state', 'S69', 'S87'], {'state': ['C', 'M'], 'S69':['U', 'P'], 'S87':['U', 'P']}) #S69: phosphorylation site for Erk1/Erk2
    Monomer('Puma', ['bf', 'state'], {'state': ['C', 'M']})

    # **Cytochrome C and Smac**
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})

def apoptosis_initial():
    alias_model_components()
    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bad(bf=None, state='C', S75='U', S99='U', S112='U'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C', S184='U'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None, state='C', S70='U'), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(Noxa(bf=None, state='C'), Noxa_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)
    Initial(Bim(bf=None, state='C', S69='U', S87='U'), Bim_0)
    Initial(Puma(bf=None, state='C'), Puma_0)
    Initial(C8(bf=None, state='pro'), C8_0)
    
def translocate_tBid_Bax_BclxL():
    """tBid, Bax and BclXL translocate to the mitochondrial membrane.
    This version of module needed instead of lopez_modules.embedded version because Bax can be phosphorylated by AKT and will not be activated."""
    alias_model_components()
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])

    free_Bax = Bax(bf=None, s1=None, s2=None, S184='U') # Alias for readability
    equilibrate(free_Bax(state='C'), free_Bax(state='M'),
                transloc_rates)

    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                transloc_rates) 
   
def apoptosis_sensitizer_translocation():
    """Translocation of apoptotic sensitizers from cytosol to mitochondria."""
    
    equilibrate(Bad(state='C', S75='U', S99='U', S112='U', bf=None), Bad(state='M', S75='U', S99='U', S112='U', bf=None), par['Bad_cyto_to_mito'])
    equilibrate(Noxa(state='C', bf=None), Noxa(state='M', bf=None), par['Noxa_cyto_to_mito'])
    equilibrate(Bim(state='C', S69='U', S87='U', bf=None), Bim(state='M', S69='U', S87='U', bf=None), par['Bim_cyto_to_mito'])
    equilibrate(Puma(state='C', bf=None), Puma(state='M', bf=None), par['Puma_cyto_to_mito'])
    
def apoptosis_bim_and_puma_bind_anti_apoptotics():
    """Binding of Bim and Puma to anti-apoptotic proteins Bcl2, Bcl-XL, and Mcl1."""
    
    bind_table([[                        Bcl2(S70='U'),                      BclxL(state='M'),           Mcl1(state='M')],
                [Bim(state='M'),        (par['Bim_bind_Bcl2']),     (par['Bim_bind_BclXL']),    (par['Bim_bind_Mcl1'])],
                [Puma(state='M'),       (par['Puma_bind_Bcl2']),    (par['Puma_bind_BclXL']),   (par['Puma_bind_Mcl1'])]],
               'bf', 'bf')

def apoptosis_bim_activate_bax():
    """Isoforms BimS and Bim-alpha3 can interact with Bax and activate it.
    Reference:
    Sarosiek, K. A. et al. BID preferentially activates BAK while BIM preferentially activates BAX, affecting chemotherapy response. Molecular Cell 51, 751–765 (2013).
    """
    
    catalyze_state(Bim(state='M'), 'bf', Bax(s1=None, s2=None), 'bf', 'state', 'M', 'A', (par['Bim_activate_Bax']))
    
def effectors_bind_anti_apoptotics():
    """
    Slightly modified from lopez_modules.embedded version to add that Bcl2 cannot be phosphorylated on S70 and have anti-apoptotic activity. (These interactions turned on in crosstalk modules).    
    Binding of Bax and Bak to Bcl2, BclxL, and Mcl1.

    Affinities of Bak for Bcl-xL and Mcl-1 are taken from Willis et al.

    Preferential affinity of Bax for Bcl-2 and Bcl-xL were taken from Zhai et
    al.  Bax:Bcl2 and Bax:Bcl-xL affinities were given order of magnitude
    estimates of 10nM.

    See comments on units for :py:func:`tBid_binds_all_anti_apoptotics`.

    Willis, S. N., Chen, L., Dewson, G., Wei, A., Naik, E., Fletcher, J. I.,
    Adams, J. M., et al. (2005). Proapoptotic Bak is sequestered by Mcl-1 and
    Bcl-xL, but not Bcl-2, until displaced by BH3-only proteins. Genes &
    Development, 19(11), 1294-1305. `doi:10.1101/gad.1304105`

    Zhai, D., Jin, C., Huang, Z., Satterthwait, A. C., & Reed, J. C. (2008).
    Differential regulation of Bax and Bak by anti-apoptotic Bcl-2 family
    proteins Bcl-B and Mcl-1. The Journal of biological chemistry, 283(15),
    9580-9586.  `doi:10.1074/jbc.M708426200`
    """

    alias_model_components()
    bind_table([[                            Bcl2(S70='U'),  BclxL(state='M'),         Mcl1],
                [Bax(active_monomer), 10e-9*N_A*V,       10e-9*N_A*V,         None],
                [Bak(active_monomer),        None,       50e-9*N_A*V,  10e-9*N_A*V]],
               kf=1e6/(N_A*V))
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 18:29:46 2015

@author: Erin
"""

from pysb import *
from pysb.macros import *
from pysb.util import alias_model_components

from parameter_dict_A431 import parameter_dict as par

def crosstalk_mapk_akt_monomers():
    Monomer('Pase9t', ['bgab1'])
    
def crosstalk_mapk_akt_initial():
    alias_model_components()
    Initial(Pase9t(bgab1=None), Pase9t_0)

def crosstalk_mapk_akt_events():
    """Defines crosstalk events between MAPK and AKT pathways.
    References:
    Kodaki, T., Woscholski, R., Hallberg, B., Rodriguez-Viciana Julian Downward, P. & Parker, P. J. The activation of phosphatidylinositol 3-kinase by Ras. Curr Biol 4, 798–806 (1994).
    Vanhaesebroeck, B., Stephens, L. & Hawkins, P. PI3K signalling: the path to discovery and understanding. Nat Rev Mol Cell Biol 13, 195–203 (2012).
    """
    
    alias_model_components()
    #ERK:P:P phosphorylates ErbBdimer-GRB2-GAB1:P (making it unable to bind PI3K)
    catalyze_state(ERK(st='PP', loc='C'), 'b', GAB1(bgrb2=ANY, bshp2=None, bpi3k=None), 'bERKPP', 'S', 'P', 'PP', (par['ERKPP_phos_GAB1P']))

    #ErbBdimer-GRB2-GAB1:P:P is dephosphorylated by Pase9t
    catalyze_state(Pase9t(), 'bgab1', GAB1(bgrb2=ANY), 'bPase9t', 'S', 'PP', 'P', (par['Pase9t_dephos_GAB1PP']))

    #AKT:PP phosphorylates RAF:P at Ser295, preventing MEK phosphorylation.
    catalyze_state(AKT(S='PP', bpip3=None), 'bpdk1', RAF(st='P'), 'b', 'ser259', 'U', 'P', (par['AKTPP_phos_RAFP']))

    #RAS-GTP binds PI3K and activates PI3K catalytic function.
    bind(RAS(bsos=None, braf=None, st='GTP'), 'bpi3k', PI3K(bgab1=None, bpip=None, berb=None), 'bras', par['RASGTP_bind_PI3K'])
    
    catalyze_state(PI3K(bras=ANY, bgab1=None, berb=None), 'bpip', PIP(bpdk1=None, bakt=None), 'bpi3k', 'S', 'PIP2', 'PIP3', par['RAS_PI3K_cat_PIP'])
    
def crosstalk_erbb_apoptosis_monomers():
    Monomer('FOXO', ['active', 'loc', 'b'], {'active':['Y', 'N'], 'loc':['C', 'N']}) #FOXO 1/3 transcription factors
    Monomer('PP2A', ['b'])

def crosstalk_erbb_apoptosis_initial():
    alias_model_components()
    Initial(FOXO(active='Y', loc='C', b=None), FOXO_0)

def crosstalk_erbb_apoptosis_events(akt_puma=True, erk_bim=True, s6k_bad=True, akt_bad=True, akt_bax=True, akt_bim=True, rsk_bad=False, erk_bcl2=True):
    """Crosstalk interactions between pathways downstream of ErbB signaling and apoptotic signaling.
    References:
    Bean, G. R. et al. PUMA and BIM Are Required for Oncogene Inactivation-Induced Apoptosis. Science Signaling 6, ra20–ra20 (2013).
    Ley, R., Ewings, K. E., Hadfield, K. & Cook, S. J. Regulatory phosphorylation of Bim: sorting out the ERK from the JNK. Cell Death & Differentiation 12, 1008–1014 (2005).
    Datta, S. R. et al. Akt phosphorylation of BAD couples survival signals to the cell-intrinsic death machinery. Cell 91, 231–241 (1997).
    Gardai, S. J. et al. Phosphorylation of Bax Ser184 by Akt Regulates Its Activity and Apoptosis in Neutrophils. Journal of Biological Chemistry 279, 21085–21095 (2004).
    Qi, X.-J., Wildey, G. M. & Howe, P. H. Evidence that Ser87 of BimEL is phosphorylated by Akt and regulates BimEL apoptotic function. J Biol Chem 281, 813–823 (2006).
    Tan, Y., Ruan, H., Demeter, M. R. & Comb, M. J. p90RSK Blocks Bad-mediated Cell Death via a Protein Kinase C-dependent Pathway. Journal of Biological Chemistry 274, 34859–34867 (1999).
    Tamura, Y., Simizu, S. & Osada, H. The phosphorylation status and anti-apoptotic activity of Bcl-2 are regulated by ERK and protein phosphatase 2A on the mitochondria. FEBS Letters 569, 249–255 (2004).    
    """
    
    alias_model_components()
    if akt_puma:
        Rule('FOXO_cyto_to_nucleus',
             FOXO(active='Y', loc='C', b=None) <>
             FOXO(active='Y', loc='N', b=None),
             *par['FOXO_cyto_to_nucleus'])
        
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', FOXO(loc='C'), 'b', 'active', 'Y', 'N', par['Akt_prevent_FOXO_transport'])
        
        #This is almost certainly not the best way to handle this but it's somewhere to start
        Rule('FOXO_tf',
             FOXO(loc='N') >>
             FOXO(loc='N') + Puma(bf=None, state='C'),
             par['FOXO_tf_Puma'])
        
    if erk_bim:
        catalyze_state(ERK(st='PP', loc='C'), 'b', Bim(state='C'), 'bf', 'S69', 'U', 'P', par['Erk_phos_Bim'])
        
        degrade(Bim(state='C', bf=None, S69='P'), par['kdeg_Bim'])
        
    if s6k_bad:
        catalyze_state(S6K(T252='P', T412='P'), 'b', Bad(state='C'), 'bf', 'S75', 'U', 'P', par['S6K_phos_Bad'])
    
    if akt_bad:
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', Bad(state='C'), 'bf', 'S99', 'U', 'P', par['Akt_phos_Bad'])
        
    if akt_bax:
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', Bax(state='C'), 'bf', 'S184', 'U', 'P', par['Akt_phos_Bax'])
    
    if akt_bim:
        catalyze_state(AKT(S='PP', bpdk1=None), 'bpip3', Bim(state='C'), 'bf', 'S87', 'U', 'P', par['Akt_phos_Bim'])
    
    if rsk_bad:
        #This is thought to be PKC dependent, so set to not run by default.
        catalyze_state(RSK1(S380='P', T573='P', S221='P'), 'b', Bad(state='C'), 'bf', 'S112', 'U', 'P', par['Rsk_phos_Bad'])
        
        degrade(Bad(state='C', bf=None, S112='P'), par['kdeg_Bad'])
    
    if erk_bcl2:
       catalyze_state(ERK(st='PP', loc='C'), 'b', Bcl2(state='C'), 'bf', 'S70', 'U', 'P', par['Erk_phos_Bcl2'])
       
       catalyze_state(PP2A(), 'b', Bcl2(state='M'), 'bf', 'S70', 'P', 'U', par['PP2A_dephos_Bcl2'])
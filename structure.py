from pymatgen.core import Structure, Element
from pymatgen.io.cif import CifWriter
import random
import numpy as np


import re
from collections import defaultdict
O_PATTERN = "O[\-\+\.\=A-Za-rt-z0-9]*$"

ALL_ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
from pymatgen.core.composition import CompositionError, Composition
def chem_dict_from_pymatgen(string):
    '''Get dictionary with chemical composition from pymatgen. Ignore weird behaviour of pymatgen when dealing with non existing elements.'''
    chem_dict = dict(Composition(string, strict=True).as_dict())
    # Fix weird behaviour of pymatgen
    for el in list(chem_dict.keys()):
        if "0+" in el:
            val = chem_dict[el]
            del chem_dict[el]
            el = el.replace("0+", "")
            chem_dict[el] = val
    return(chem_dict)
import re
def isfloat(value):
    """"Check if variable is float."""
    try:
      float(value)
      return True
    except ValueError:
      return False
def get_chem_dict(string, weird_oxygen=False):
    """Get a sorted dictionary of the elements in the string and the corresponding values.
    
    Can deal with unspecified oxygen quantitites if they follow the pattern O_PATTERN. If oxygen is not specified at all (e.g. Oz) it's value is np.nan, if it's value is something like O7+z it's new value is 7. The sorting is alphabetically for the keys, except that oxygen is always last.
    """
    # Stupid exception in case someone used Y as unknown value for oxygen.
    exception_pattern = "OY$"
    if re.search(exception_pattern, string) != None:
        O_value = "Y"
        string = re.sub(exception_pattern, "", string)
    
    try:
        chem_dict = chem_dict_from_pymatgen(string)
    except (CompositionError, ValueError) as error: 
        # Search for weird oxygen strings in the Supercon.
        if weird_oxygen == True and re.search(O_PATTERN, string) != None:
            split_idx = re.search(O_PATTERN, string).start()
            O_string = string[split_idx:]
            O_value = O_string.split("O")[-1]
            if O_value == "nan":
                O_value = np.nan
            string = string[:split_idx]
        elif "Onan" in string:
            # If value of oxygen is np.nan.
            O_value = np.nan
            string = string.replace("Onan", "")
        try:
            chem_dict = chem_dict_from_pymatgen(string)
        except (CompositionError, ValueError) as error:
            # print("Bad formula, returned empty chemdict: {}".format(string))
            chem_dict = {}
            return(chem_dict)
    
    # Get oxygen value if no weird string is present.
    try:
        O_value
    except NameError:
        if "O" in chem_dict.keys():
            O_value = chem_dict["O"]
            del chem_dict["O"]
        else:
            O_value = None

    assert all([isfloat(chem_dict[key]) for key in chem_dict.keys()]), "In the chem_dict we found an element with not float value: {}".format(chem_dict)
    
    # Sort alphabetically, except oxygen comes last.
    sorted_elements = sorted(chem_dict.keys())
    sorted_chem_dict = {el: chem_dict[el] for el in sorted_elements}
    if O_value != None:
        sorted_chem_dict["O"] = O_value
    return(sorted_chem_dict)


def normalise_quantities(chemdict, total=100):
    """Normalises the quantities in the chemical formula dictionary so that the new values are in percent of the sum of quantities. Can deal with np.nan as quantity of an element.
    """
    if "O" in chemdict.keys():
        # Oxygen value may either be number or np.nan
        assert isfloat(chemdict["O"]), chemdict
    quantities = np.array(list(chemdict.values()))
    total_quantities = np.nansum(quantities)
    norm_factor = total/total_quantities
    chemdict = {el: norm_factor*chemdict[el] for el in chemdict}
    return(chemdict)
def standardise_chem_formula(string, normalise=False):
    '''Does some normalising and sorting, so that chemical formulas have the same string when they are chemically identical.'''
    string = string.replace(" ", "")
    chem_dict = get_chem_dict(string)
    if normalise == True:
        chem_dict = normalise_quantities(chem_dict)
            
    if "O" in chem_dict.keys():
        O_value = chem_dict["O"]
        del chem_dict["O"]
    else:
        O_value = ""
    
    chem_list = [el + str(round(chem_dict[el], 3)) for el in chem_dict.keys()]
    if O_value != "":
        chem_list.append("O{}".format(O_value))
    chem_list = [re.sub("\.0$", "", el) for el in chem_list]
    sorted_string = "".join(chem_list)
    return(sorted_string)

def parse_formula(formula):
    # 匹配化学式中的元素和数量
    pattern = re.compile(r'([A-Z][a-z]*)(\d*\.?\d*)')
    matches = pattern.findall(formula)
    
    # 解析并计算每种元素的数量
    elements = defaultdict(float)
    for element, count in matches:
        if count == '':
            count = 1  # 如果没有给出数量，则默认为1
        else:
            count = float(count)
        elements[element] += count
    
    return dict(elements)

def compare_formulas(formula1, formula2):
    # 解析两个化学式
    elements1 = parse_formula(formula1)
    elements2 = parse_formula(formula2)
    
    # 获取所有的元素集合
    all_elements = set(elements1.keys()).union(set(elements2.keys()))
    
    # 打印对比结果
    comparison = []
    for element in all_elements:
        count1 = elements1.get(element, 0)
        count2 = elements2.get(element, 0)
        comparison.append((element, count1, count2))
    
    return comparison
def compare_structures(target_formula,reference_formula,hold=0.189):
    comparison = compare_formulas(target_formula, reference_formula)

    div_min=min([abs(count2-count1) for (element, count1, count2) in comparison])
    comparison_tmp=comparison
    if not div_min<0.07:
        comparison_tmp=[(element, count1*2, count2) for (element, count1, count2) in comparison]
        div_min=min([abs(count2-count1) for (element, count1, count2) in comparison_tmp])
        if not div_min<0.07:
            comparison_tmp=[(element, count1*3, count2) for (element, count1, count2) in comparison]
            div_min=min([abs(count2-count1) for (element, count1, count2) in comparison_tmp])
            if not div_min<0.07:
                comparison_tmp=[(element, count1*4, count2) for (element, count1, count2) in comparison]
                div_min=min([abs(count2-count1) for (element, count1, count2) in comparison_tmp])
    comparison=comparison_tmp
    div_max=max([abs(count2-count1) for (element, count1, count2) in comparison])
    if div_max<hold:
        return True
    else:
        return False

#compare_structures('Ba2Cu2.875Y1Zn0.125O6.93','Ba2Cu2.875Y1Zn0.125O6.93') 
from pymatgen.core import Composition
def get_doping_structure(target_formula,structure):
    reference_formula=structure.formula.replace(" ","")
    comparison = compare_formulas(target_formula, reference_formula)
    #print(comparison)
    sumcount1=sum([count1 for (element, count1, count2) in comparison])
    sumcount2=sum([count2 for (element, count1, count2) in comparison])

    comparison=[(element, count1*round(sumcount2/sumcount1), count2) for (element, count1, count2) in comparison]

    max_div=[count2-count1 for (element, count1, count2) in comparison]
    max_ele=max_div.index(max(max_div))
    element, count1, count2=comparison[max_ele]
    supercell=structure
    if (count2-count1)>0.75:
        supercell=structure
    elif (count2-count1)>0.45:    
        supercell = structure * (1, 1, 2)
        comparison=[(element, count1*2, count2*2) for (element, count1, count2) in comparison]
    elif(count2-count1)>0.29:
        supercell = structure * (1, 1, 3)
        comparison=[(element, count1*3, count2*3) for (element, count1, count2) in comparison]
    elif(count2-count1)>0.19:
        supercell = structure * (2, 2, 1)
        comparison=[(element, count1*4, count2*4) for (element, count1, count2) in comparison]
    elif(count2-count1)>0.10:
        supercell = structure * (2, 2, 2)
        comparison=[(element, count1*8, count2*8) for (element, count1, count2) in comparison]
    else:
        return structure

    elements=[element for (element, count1, count2) in comparison]
    if "O" in elements and "F" not in elements:
        element, count1, count2=comparison[elements.index('O')]
        del comparison[elements.index('O')]
        if count2-count1>0.75:
            num_to_replace=round(count2-count1)
            r_indices = [i for i, site in enumerate(supercell.sites) if Composition(site.species).chemical_system == 'O']
            #i=random.choice(r_indices)
            indices_to_replace = random.sample(r_indices, num_to_replace)
            supercell.remove_sites(indices_to_replace)

    max_div=np.around([count2-count1 for (element, count1, count2) in comparison],3).tolist()
    if max(max_div)<0.19:
        return supercell
    max_ele=max_div.index(max(max_div))
    element1, count1, count2=comparison[max_ele]
    min_ele=max_div.index(min(max_div))
    element2, _, _=comparison[min_ele]

    r_indices = [i for i, site in enumerate(supercell.sites) if Composition(site.species).chemical_system == element1]
    if len(r_indices)<1:
        return supercell
    num_to_replace=round(count2-count1)
    indices_to_replace = random.sample(r_indices, num_to_replace)  
    for i in indices_to_replace:
        supercell.replace(i, element2)
    #i=random.choice(r_indices)
    #supercell.replace(i, element2)
    return supercell


import copy

import numpy as np
import pandas as pd
def similarity_chem_formula(chem_dict_sc, chem_dict_2, max_relcutoff, total_relcutoff, min_abscutoff):
    """Checks if the chemical formula is "similar" as defined. Assumes that chemdicts have the same elements and are sorted.
    """
    # To not modify the actual dictionary from out of this function.
    chemdict_sc, chemdict_2 = copy.deepcopy(chem_dict_sc), copy.deepcopy(chem_dict_2)
    
    assert chemdict_sc.keys() == chemdict_2.keys()
        
    quantities_sc = np.array(list(chemdict_sc.values()))
    quantities_2 = np.array(list(chemdict_2.values()))
    # Normalise in case one formula is conventional unit cell and the other formula is primitive unit cell.
    norm = np.sum(quantities_sc)/np.sum(quantities_2)
    quantities_2 = quantities_2*norm
    
    # Check differences of elements.
    diffs, reldiffs, totreldiff = calculate_numeric_deviation(quantities_sc, quantities_2)
    
    # First condition
    if totreldiff > total_relcutoff:
        return(False, totreldiff)
    
    # Second condition
    for diff, reldiff in zip(diffs, reldiffs):
        if diff > min_abscutoff and reldiff > max_relcutoff:
            return(False, totreldiff)
        
    return(True, totreldiff)


def calculate_numeric_deviation(nums1, nums2):
    """Calculates absolute and relative differences.
    """
    
    # Make sure nums are numpy arrays.
    nums1, nums2 = np.array(nums1), np.array(nums2)
    
    # Calculate differences.
    diffs = np.abs(nums1 - nums2)
    # Calculate relative differences.
    reldiffs = 2*diffs/(nums1 + nums2)
    # The formula for totreldiff is not the sum of reldiffs but rather the weighted sum of reldiffs. That means elements with low quantities are weighed less than elements with high quantities. This is the right behaviour, because otherwise small changes of doping of small amounts would dominate (e.g. if element1 has quantity 0.1 and element2 has quantity 0.15 they would already make a difference of about 40%, this should be weighed down by other elements with higher quantities.)
    totreldiff = 2*diffs.sum()/(nums1.sum() + nums2.sum())
    return(diffs, reldiffs, totreldiff)

n_max_doping_elements = 1

# For similarity == 2
lower_max_relcutoff = 0.10001
lower_total_relcutoff = 0.05001
lower_min_abscutoff = 0.15001

# For similarity == 3
higher_max_relcutoff = 0.20001
higher_total_relcutoff = 0.15001
higher_min_abscutoff = 0.3001

def get_formula_similarity(chemdict_sc, chemdict_2, 
                           lower_max_relcutoff=0.10001, 
                           lower_total_relcutoff= 0.05001, 
                           lower_min_abscutoff= 0.15001, 
                           higher_max_relcutoff= 0.20001, 
                           higher_total_relcutoff= 0.15001, 
                           higher_min_abscutoff= 0.3001):
    """Get a score for the similarity of the chemical formulas.
    """    
    elements_sc = list(chemdict_sc.keys())
    elements_2 = list(chemdict_2.keys())
    
    # Add additonally doped elements from the Supercon formula into the cif formula (i.e. pad chemical formula with zeros).
    assert len(elements_sc) >= len(elements_2)
    assert all([el in elements_sc for el in elements_2])
    #if len(elements_sc) > len(elements_2):
    #    if len(elements_2) == 1:
            # If the cif formula is an element we don't want any doped matches because it's likely that even a bit doping changes the structure a lot.
    #        formula_similarity = np.nan
    #        totreldiff = np.nan
    #        return(formula_similarity, totreldiff)
    chemdict_2 = {el: chemdict_2[el] if el in elements_2 else 0. for el in elements_sc}
    
    # Check if chemical formulas are very close ('similar')
    formulas_similar, totreldiff = similarity_chem_formula(
                                        chem_dict_sc = chemdict_sc,
                                        chem_dict_2 = chemdict_2,
                                        max_relcutoff = lower_max_relcutoff,
                                        total_relcutoff = lower_total_relcutoff,
                                        min_abscutoff = lower_min_abscutoff
                                        )

    if totreldiff == 0:
        # relative formulas are the same.
        formula_similarity = 1
    elif formulas_similar == True:
        formula_similarity = 2
    else:
        formula_similarity = 4

    # Check if chemical formulas are relatively close ('doped').
    if pd.isna(formula_similarity): 
        formulas_doped, totreldiff = similarity_chem_formula(
                                        chem_dict_sc = chemdict_sc,
                                        chem_dict_2 = chemdict_2,
                                        max_relcutoff = higher_max_relcutoff,
                                        total_relcutoff = higher_total_relcutoff,
                                        min_abscutoff = higher_min_abscutoff
                                        )
        if formulas_doped == True:
            formula_similarity = 3
        
    return(formula_similarity, totreldiff)

def doping(structure,min_dop=0.4):
    i=0
    for site in structure:
        # site.species_and_occu 是一个字典，键是元素，值是该元素的占位（掺杂浓度）
        #print(site.species.items())
        for element, occupancy in site.species.items():
            if occupancy > min_dop and occupancy<0.5:
                i=i+1
    if i==1:
        return True
    else:
        return False
from pymatgen.core.sites import PeriodicSite
def doping_occupancy(structure,min_dop=0.25):
    equivalent_sites: list[list[int]] = []
    exemplars: list[PeriodicSite] = []
    for idx, site in enumerate(structure):
        if site.is_ordered:
            continue
        for j, ex in enumerate(exemplars):
            sp = ex.species
            if not site.species.almost_equals(sp):
                continue
            equivalent_sites[j].append(idx)
            break
        else:
            equivalent_sites.append([idx])
            exemplars.append(site)
    for group in equivalent_sites:
        total_occupancy = dict(
            sum((structure[idx].species for idx in group), Composition()).items()  # type: ignore[attr-defined]
        )
        # round total occupancy to possible values
        for key, val in total_occupancy.items():
            if abs(val - round(val)) > min_dop:
                return True
    return False

def doping_occupancys(structure,min_dop=0.25):
    equivalent_sites: list[list[int]] = []
    exemplars: list[PeriodicSite] = []
    for idx, site in enumerate(structure):
        if site.is_ordered:
            continue
        for j, ex in enumerate(exemplars):
            sp = ex.species
            if not site.species.almost_equals(sp):
                continue
            equivalent_sites[j].append(idx)
            break
        else:
            equivalent_sites.append([idx])
            exemplars.append(site)

    doping_occupancys=0
    for group in equivalent_sites:
        total_occupancy = dict(
            sum((structure[idx].species for idx in group), Composition()).items()  # type: ignore[attr-defined]
        )
        # round total occupancy to possible values
        for key, val in total_occupancy.items():
            if abs(val - round(val)) > doping_occupancys:
                doping_occupancys=abs(val - round(val))
    return doping_occupancys

def check_for_doping(structure):

    for i, site in enumerate(structure):
        if len(site.species) > 1:
            return True
    return False
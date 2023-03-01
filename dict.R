### UMLS-HPO ANALYSIS OF INDIVIDUALS WITH GENETIC EPILEPSY SYNDROMES -----------
## Lookup tables
##
## Author: Christian Bosselmann, MD
##
## Date Created: 2023-02-07
##
## Copyright (c) Christian Bosselmann, 2023
## Email: bosselc@ccf.org
##
### ----------------------------------------------------------------------------

### manual recoding of ASMs 
recodeASM <- function(df_asm){
  vec_asm <- vector()
  vec_asm <- df_asm %>%
    mutate(MED_NAME = case_when(str_detect(MED_NAME, regex('vimpat', ignore_case = T)) ~ 'lacosamide',
                                str_detect(MED_NAME, regex('lacosamide', ignore_case = T)) ~ 'lacosamide', 
                                str_detect(MED_NAME, regex('banzel', ignore_case = T)) ~ 'rufinamide',
                                str_detect(MED_NAME, regex('tegretol', ignore_case = T)) ~ 'carbamazepine',
                                str_detect(MED_NAME, regex('carbatrol', ignore_case = T)) ~ 'carbamazepine',
                                str_detect(MED_NAME, regex('topamax', ignore_case = T)) ~ 'topiramate',
                                str_detect(MED_NAME, regex('topiram', ignore_case = T)) ~ 'topiramate',
                                str_detect(MED_NAME, regex('qudexy', ignore_case = T)) ~ 'topiramate',
                                str_detect(MED_NAME, regex('trokendi', ignore_case = T)) ~ 'topiramate',
                                str_detect(MED_NAME, regex('eprontia', ignore_case = T)) ~ 'topiramate',
                                str_detect(MED_NAME, regex('spritam', ignore_case = T)) ~ 'levetiracetam',
                                str_detect(MED_NAME, regex('keppra', ignore_case = T)) ~ 'levetiracetam',
                                str_detect(MED_NAME, regex('levetir', ignore_case = T)) ~ 'levetiracetam',
                                str_detect(MED_NAME, regex('sympazan', ignore_case = T)) ~ 'clobazam',
                                str_detect(MED_NAME, regex('onfi', ignore_case = T)) ~ 'clobazam',
                                str_detect(MED_NAME, regex('loraze', ignore_case = T)) ~ 'lorazepam',
                                str_detect(MED_NAME, regex('midazol', ignore_case = T)) ~ 'midazolam',
                                str_detect(MED_NAME, regex('nayzilam', ignore_case = T)) ~ 'midazolam',
                                str_detect(MED_NAME, regex('klonopin', ignore_case = T)) ~ 'clonazepam',
                                str_detect(MED_NAME, regex('lyrica', ignore_case = T)) ~ 'pregabalin',
                                str_detect(MED_NAME, regex('pregabalin', ignore_case = T)) ~ 'pregabalin',
                                str_detect(MED_NAME, regex('zonis', ignore_case = T)) ~ 'zonisamide',
                                str_detect(MED_NAME, regex('zoneg', ignore_case = T)) ~ 'zonisamide',
                                str_detect(MED_NAME, regex('valtoco', ignore_case = T)) ~ 'diazepam',
                                str_detect(MED_NAME, regex('diastat', ignore_case = T)) ~ 'diazepam',
                                str_detect(MED_NAME, regex('valpro', ignore_case = T)) ~ 'valproic acid',
                                str_detect(MED_NAME, regex('depak', ignore_case = T)) ~ 'valproic acid',
                                str_detect(MED_NAME, regex('lacosamide', ignore_case = T)) ~ 'lacosamide',
                                str_detect(MED_NAME, regex('fycompa', ignore_case = T)) ~ 'perampanel',
                                str_detect(MED_NAME, regex('perampanel', ignore_case = T)) ~ 'perampanel',
                                str_detect(MED_NAME, regex('fintepla', ignore_case = T)) ~ 'fintepla',
                                str_detect(MED_NAME, regex('fenfluramine', ignore_case = T)) ~ 'fenfluramine',
                                str_detect(MED_NAME, regex('lamot', ignore_case = T)) ~ 'lamotrigine',
                                str_detect(MED_NAME, regex('oxc', ignore_case = T)) ~ 'oxcarbazepine',
                                str_detect(MED_NAME, regex('trileptal', ignore_case = T)) ~ 'oxcarbazepine',
                                str_detect(MED_NAME, regex('loraz', ignore_case = T)) ~ 'lorazepam',
                                str_detect(MED_NAME, regex('vigabatr', ignore_case = T)) ~ 'vigabatrin',
                                str_detect(MED_NAME, regex('sabril', ignore_case = T)) ~ 'vigabatrin',
                                str_detect(MED_NAME, regex('tiagab', ignore_case = T)) ~ 'tiagabine',
                                str_detect(MED_NAME, regex('rufinam', ignore_case = T)) ~ 'rufinamide',
                                str_detect(MED_NAME, regex('primid', ignore_case = T)) ~ 'primidone',
                                str_detect(MED_NAME, regex('phenyt', ignore_case = T)) ~ 'phenytoin',
                                str_detect(MED_NAME, regex('dilantin', ignore_case = T)) ~ 'phenytoin',
                                str_detect(MED_NAME, regex('epanutin', ignore_case = T)) ~ 'phenytoin',
                                str_detect(MED_NAME, regex('cerebyx', ignore_case = T)) ~ 'phenytoin',
                                str_detect(MED_NAME, regex('phenobarb', ignore_case = T)) ~ 'phenobarbital',
                                str_detect(MED_NAME, regex('gabapent', ignore_case = T)) ~ 'gabapentin',
                                str_detect(MED_NAME, regex('felbamat', ignore_case = T)) ~ 'felbamate',
                                str_detect(MED_NAME, regex('everol', ignore_case = T)) ~ 'everolimus',
                                str_detect(MED_NAME, regex('ethosux', ignore_case = T)) ~ 'ethosumixide',
                                str_detect(MED_NAME, regex('zaront', ignore_case = T)) ~ 'ethosumixide',
                                str_detect(MED_NAME, regex('cenoba', ignore_case = T)) ~ 'cenobamate',
                                str_detect(MED_NAME, regex('xcopri', ignore_case = T)) ~ 'cenobamate',
                                str_detect(MED_NAME, regex('cannab', ignore_case = T)) ~ 'cannabidiol',
                                str_detect(MED_NAME, regex('epidiolex', ignore_case = T)) ~ 'cannabidiol',
                                str_detect(MED_NAME, regex('briv', ignore_case = T)) ~ 'brivaracetam',
                                str_detect(MED_NAME, regex('stiri', ignore_case = T)) ~ 'stiripentol',
                                str_detect(MED_NAME, regex('celontin', ignore_case = T)) ~ 'mesuximide',
                                str_detect(MED_NAME, regex('eslicarbazepine', ignore_case = T)) ~ 'eslicarbazepine acetate')) %>%
    pull(MED_NAME)
  return(vec_asm)
}

### code ASMs by mechanism subgroup
# define ASM groups, cf. doi.org/10.1007/s40263-021-00827-8
groupASM <- function(asm_vec){
  asm_map <- tibble()
  asm_map <- tibble(MED_NAME = asm_vec,
                    MED_GROUP = NA) %>%
    mutate(MED_GROUP = case_when(
      str_detect(MED_NAME, regex('felbam', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('valpr', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('hormon', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('zonis', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('rufin', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('cannab', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('cenob', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('topira', ignore_case = T)) ~ 'Mixed/unknown',
      str_detect(MED_NAME, regex('phenyt', ignore_case = T)) ~ 'Na',
      str_detect(MED_NAME, regex('carbamaz', ignore_case = T)) ~ 'Na',
      str_detect(MED_NAME, regex('carbaz', ignore_case = T)) ~ 'Na',
      str_detect(MED_NAME, regex('lamotr', ignore_case = T)) ~ 'Na',
      str_detect(MED_NAME, regex('lacosam', ignore_case = T)) ~ 'Na',
      str_detect(MED_NAME, regex('suxim', ignore_case = T)) ~ 'Ca',
      str_detect(MED_NAME, regex('gabap', ignore_case = T)) ~ 'Ca',
      str_detect(MED_NAME, regex('pregabal', ignore_case = T)) ~ 'Ca',
      str_detect(MED_NAME, regex('phenobarb', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('primidon', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('stirip', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('azol', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('azepam', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('clobaz', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('tiaga', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('vigabatr', ignore_case = T)) ~ 'GABA',
      str_detect(MED_NAME, regex('peramp', ignore_case = T)) ~ 'AMPA',
      str_detect(MED_NAME, regex('brivara', ignore_case = T)) ~ 'SV2A',
      str_detect(MED_NAME, regex('levetir', ignore_case = T)) ~ 'SV2A',
      str_detect(MED_NAME, regex('fenfl', ignore_case = T)) ~ '5-HT',
      str_detect(MED_NAME, regex('everol', ignore_case = T)) ~ 'MTOR'
    ) )
  return(asm_map)
}

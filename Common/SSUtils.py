
## the below 6 lines shall be moved into another script. Need to double check conversion from SS8 to SS3
SS8Letter2Code = {'H':0, 'G':1, 'I':2, 'E':3, 'B':4, 'T':5, 'S':6, 'L':7, 'C':7 }
SS8Letter2SS3Code = {'H':0, 'G':0, 'I':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2, 'C':2 }
SS8Letter2SS3Letter = {'H':'H', 'G':'H', 'I':'H', 'E':'E', 'B':'E', 'T':'L', 'S':'L', 'L':'L', 'C':'L' }
## in the tgt file, the predicted 3-state secondary structure is arranged in the order of H, E, and L
## in the tgt file, the predicted 8-state secondary structure is arranged in the order of H     G     I     E     B     T     S     L
## in the tgt file, the predicted solvent accessibility is arranged in the order of Bury   Medium  Exposed


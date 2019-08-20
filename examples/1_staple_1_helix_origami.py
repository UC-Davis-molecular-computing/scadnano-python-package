import scadnano as sc

def main():
    length = 9
    helices = [sc.Helix(0, length)]
    stap_ss = sc.Substrand(0, sc.right, 0, length)
    scaf_ss = sc.Substrand(0, sc.left, 0, length)
    stap = sc.Strand([stap_ss])
    scaf = sc.Strand([scaf_ss], color=sc.default_scaffold_color)
    strands = [stap, scaf]
    design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)
    design.add_deletion(helix_idx=0, offset=2)
    design.assign_dna(scaf, 'AACT' * (length // 4))

    return design

if not sc.in_browser() and __name__ == '__main__':
    design = main()
    design.write_file(directory='output_designs')

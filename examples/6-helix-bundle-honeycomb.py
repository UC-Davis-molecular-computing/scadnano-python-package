import scadnano as sc


def main():
    max_offset = 1295
    helices = [
        sc.Helix(grid_position=(1, 0, 0), max_offset=max_offset),
        sc.Helix(grid_position=(0, 1, 0), max_offset=max_offset),
        sc.Helix(grid_position=(1, 2, 0), max_offset=max_offset),
        sc.Helix(grid_position=(2, 2, 0), max_offset=max_offset),
        sc.Helix(grid_position=(2, 1, 0), max_offset=max_offset),
        sc.Helix(grid_position=(2, 0, 0), max_offset=max_offset),
    ]
    stap_ss = sc.Substrand(0, sc.forward, 0, 100)
    scaf_ss = sc.Substrand(0, sc.reverse, 0, 100)
    stap = sc.Strand([stap_ss])
    scaf = sc.Strand([scaf_ss], color=sc.default_scaffold_color)
    strands = [stap, scaf]
    design = sc.DNADesign(helices=helices, strands=strands, grid=sc.honeycomb)
    design.add_deletion(helix=0, offset=2)
    design.assign_dna(scaf, sc.m13_sequence)

    return design


if not sc.in_browser() and __name__ == '__main__':
    design = main()
    design.write_scadnano_file(directory='output_designs')

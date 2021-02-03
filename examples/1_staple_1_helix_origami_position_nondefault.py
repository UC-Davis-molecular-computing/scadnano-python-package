import scadnano as sc


def create_design() -> sc.Design:
    length = 9
    pitch90_group = sc.HelixGroup(pitch=90)
    helices = [
        sc.Helix(max_offset=length, major_ticks=[2, 5], position=sc.Position3D(x=1, y=2, z=3), group='pitch90')]
    stap_ss = sc.Domain(0, True, 0, length)
    scaf_ss = sc.Domain(0, False, 0, length)
    stap = sc.Strand([stap_ss])
    scaf = sc.Strand([scaf_ss], color=sc.default_scaffold_color)
    strands = [stap, scaf]
    design = sc.Design(helices=helices, groups={
                       'pitch90': pitch90_group}, strands=strands)
    design.assign_dna(scaf, 'AACT' * (length // 4))

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
